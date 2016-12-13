#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "MTGenerator.h"

int N,gaps,activeN,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,multimerN,countCollidingPairs,ODFLength,OCFMode,skipFirstIteration,saveConfigurations,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    onlyMath[2]={0,0},initMode;
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,intervalOrientations,savedConfigurationsInt,generatorStartPoint=0;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg[2],deltaR,deltaPhi,deltaV=0.1, //*sigma(=1)
       multimerS,multimerD,randomStartStep[2],detBoxMatrix,
       neighRadius,neighRadius2,neighRadiusMod,neighSafeDistance,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][3],pi=3.1415926535898;
double L,C,ROkreguOpisanego,absoluteMinimum,absoluteMinimum2,minDistance,maxDistance,VcpPerParticle;
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferMN[20],bufferMS[20],bufferMD[20],bufferFolderIndex[5],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     orientationsFileName[200]="Orientations",
     orientatCorrelFunFileName[200]="OrientatCorrelFun",
     orientationsResultsFileName[200]="OrientatRes",
     configurationsListFileName[200]="ConfigurationsList.txt",
     loadConfigurationsFileName[200]="Configurations",
     loadedJOBID[50]="j-none";
/* //matrix 1/16
int boxCellsN[200][200],boxCellsIndex[200][200][10];
double bCSize[2];
*/
/////////////////  PARTICLE functions{
typedef struct particle {
    double r[2], normR[2];  //x,y
    double phi;   //kąt mierzony od kierunku x
    int neighbours[20], neighCounter;
    //int cell[2]; //matrix 2/16
} particle;

/*void assignParticlesToCells (particle *particles, double boxMatrix[2][2]) { //matrix 3/16
    bCSize[0]=boxMatrix[0][0]/maxDistance; bCSize[1]=boxMatrix[1][1]/maxDistance;
    for (int i=0;i<bCSize[0];i++) for (int j=0;j<bCSize[1];j++) {
        boxCellsN[i][j]=0;
        for (int k=0;k<10;k++) boxCellsIndex[i][j][k]=-1;
    }
    for (int i=0;i<N;i++) {
        int indexes[2] = {(int)(particles[i].normR[0]*bCSize[0]),(int)(particles[i].normR[1]*bCSize[1])};
        boxCellsIndex[indexes[0]][indexes[1]][boxCellsN[indexes[0]][indexes[1]]++]=i;
        particles[i].cell[0]=indexes[0]; particles[i].cell[1]=indexes[1];
    }
}*/

void updateNeighbourList (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
               normalizedRY=particles[i].normR[1]-particles[j].normR[1],
               rx=particles[i].r[0]-particles[j].r[0],
               ry=particles[i].r[1]-particles[j].r[1];
        rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
        ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
        double r2=rx*rx+ry*ry;
        if (r2<neighRadius2) {
            particles[i].neighbours[particles[i].neighCounter++]=j;
            particles[j].neighbours[particles[j].neighCounter++]=i;
        }
    }
}

double normalizeAngle (double phi) {//utrzymywanie phi w zakresie (-30st,30st) = (-Pi/n,Pi/n)
    int change;
    do {
        change=0;
        if (phi<-C) {
            phi+=2*C;
            change=1;
        } else if (phi>C) {
            phi-=2*C;
            change=1;
        }
    } while (change==1);
    return phi;
}

double get2DiscsDistance (int aNum, int bNum, double dr, double aAngle, double bAngle) {
    double discAngle = aAngle+aNum*2*C;
    double aNumDiscX = cos(discAngle)*ROkreguOpisanego,
           aNumDiscY = sin(discAngle)*ROkreguOpisanego;
    discAngle = bAngle+bNum*2*C;
    double bNumDiscX = dr+cos(discAngle)*ROkreguOpisanego,
           bNumDiscY = sin(discAngle)*ROkreguOpisanego;
    double xDistance = aNumDiscX-bNumDiscX, yDistance = aNumDiscY-bNumDiscY;
    return sqrt(xDistance*xDistance+yDistance*yDistance);
}

int checkOverlaps (double dr, double aAngle, double bAngle) {
    int overlap=0;
    double startAAngle,startBAngle;
    int a,b;
    if (aAngle>absoluteMinimum2)
        if (aAngle>bAngle) {
            startAAngle=aAngle-2*C;
            a=2;
        } else {
            startAAngle=aAngle;
            a=1;
        }
    else if (aAngle>-absoluteMinimum2) {
        startAAngle=aAngle;
        a=1;
    } else
        if (aAngle<bAngle) {
            startAAngle=aAngle;
            a=2;
        } else {
            startAAngle=aAngle;
            a=1;
        }
    if (bAngle>absoluteMinimum2)
        if (a==2) {
            startBAngle=bAngle+(double)multimerN*C;
            b=1;
        } else {
            startBAngle=bAngle+((double)multimerN-2.0)*C;
            b=2;
        }
    else if (bAngle>-absoluteMinimum2) {
        startBAngle=bAngle+(double)multimerN*C;
        b=1;
    } else
        if (a==2) {
            startBAngle=bAngle+(double)multimerN*C;
            b=1;
        } else {
            startBAngle=bAngle+(double)multimerN*C;
            b=2;
        }
    for (int j=0;j<a;j++) for (int k=0;k<b;k++)
        if (get2DiscsDistance(j,k,dr,startAAngle,startBAngle)<multimerD) {
            overlap=1;
            j=a; k=b;
            break;
        }
    return overlap;
}

double minimalDistanceAnalyticalMethodHCM (int indeksPowierzchni, double aAngle, double bAngle, double a[4], double b[4]) {
    double angleA=aAngle+a[indeksPowierzchni],angleB=bAngle+b[indeksPowierzchni],
           buffer=sin(angleA)+sin(angleB);
    return (ROkreguOpisanego*(cos(angleA)+cos(angleB))+sqrt(multimerD*multimerD-ROkreguOpisanego*ROkreguOpisanego*buffer*buffer));
}

double getMinimalDistanceAnalyticalMethodForEvenHCM (double aAngle, double bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-C,-C,C,C};
    if (aAngle<-absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13=-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<0) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve14=aAngle,
               intersectionCurve34=-C-ARCSIN;
        if (bAngle>intersectionCurve14) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve34) return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve14=aAngle,
               intersectionCurve12=C-ARCSIN;
        if (bAngle>intersectionCurve12) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve14) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    } else {
        double ARCSIN=asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24=-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    }
}

double getMinimalDistanceAnalyticalMethodForOddHCM (double aAngle, double bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-2*C,0,0,2*C};
    if (aAngle<-absoluteMinimum) {
        double ARCSIN = asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13 = C-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<0) {
        double ARCSIN = asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve12 = aAngle+C,
               intersectionCurve23 = -ARCSIN;
        if (bAngle>intersectionCurve12) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve23) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN = asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve23 = -ARCSIN,
               intersectionCurve34 = aAngle-C;
        if (bAngle>intersectionCurve23) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve34) return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    } else {
        double ARCSIN = asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24 = -C-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    }
}

int checkOverlapsAnalyticalMethodHCM (double dr, double aAngle, double bAngle) {
    int overlap=0;

    aAngle=normalizeAngle(aAngle+C);
    bAngle=normalizeAngle(bAngle+C);

    if (multimerN%2==0) {
        if (dr<getMinimalDistanceAnalyticalMethodForEvenHCM(aAngle,bAngle)) overlap=1;
    } else if (dr<getMinimalDistanceAnalyticalMethodForOddHCM(aAngle,bAngle)) overlap=1;
    return overlap;
}

int createRandomGaps (particle *particles, double boxMatrix[2][2]) {
    printf("Creating %d random gaps... ",gaps);
    InitRandomMT();
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTGenerate(randomStartStep)%1000000/1000000.0*N);
            for (int j=0;j<i;j++) {
                double normalizedRX=particles[gapsIndexes[i]].normR[0]-particles[gapsIndexes[j]].normR[0],
                       normalizedRY=particles[gapsIndexes[i]].normR[1]-particles[gapsIndexes[j]].normR[1],
                       rx=particles[gapsIndexes[i]].r[0]-particles[gapsIndexes[j]].r[0],
                       ry=particles[gapsIndexes[i]].r[1]-particles[gapsIndexes[j]].r[1];
                rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                double r2=rx*rx+ry*ry, dr=sqrt(r2);
                if (dr<neighRadius) {
                    i--;
                    innerAttempt++;
                    break;
                }
                if (j+1==i) innerAttempt=0;
            }
            if (innerAttempt>100000) break;
            if (innerAttempt==0 && i+1==gaps) allReady=1;
        }
    } while (!allReady && attempt<10000000);

    if (attempt>=10000000) {
        printf("ERROR: Couldn't create %d random gaps in %d steps.\n",gaps,attempt);
        return 0;
    } else {
        bool change; do {
            change=false;
            for (int i=gaps-1;i>0;i--) {
                if (gapsIndexes[i-1]>gapsIndexes[i]) {
                    int buffer=gapsIndexes[i];
                    gapsIndexes[i]=gapsIndexes[i-1]; gapsIndexes[i-1]=buffer;
                    change=true;
                }
            }
        } while (change);
        int actualGapIndex=0;
        for (int i=0;i<activeN;i++) {
            if (actualGapIndex<gaps) while (i+actualGapIndex==gapsIndexes[actualGapIndex]) {
                actualGapIndex++;
                if (actualGapIndex==gaps) break;
            }
            for (int j=0;j<2;j++) {
                particles[i].r[j]=particles[i+actualGapIndex].r[j];
                particles[i].normR[j]=particles[i+actualGapIndex].normR[j];
            }
        }
        printf("done\n");
        return 1;
    }
}

int initPositions (particle *particles, double boxMatrix[2][2], double matrixOfParticlesSize[2], int n[2], double matrixCellXY[6][6][2], double matrixCellPhi[6][6], double pacFrac) {
    double mod=sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]), interval[2][2], actualPosition[2];
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) interval[i][j]=boxMatrix[i][j]/matrixOfParticlesSize[j]/mod*n[j];
    int rowCounter=0, columnCounter=0;
    for (int i=0;i<N;i++) {
        int cellNumber[2][2]={{columnCounter/n[0],rowCounter/n[1]},{columnCounter%n[0],rowCounter%n[1]}}; //cellNumber[0/1][X/Y]: 0-numer komorki, 1-kolumna/rzad W komorce
        actualPosition[0]=cellNumber[0][0]*interval[0][0]+cellNumber[0][1]*interval[0][1]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][0]*sqrt(pacFrac);
        actualPosition[1]=cellNumber[0][1]*interval[1][1]+cellNumber[0][0]*interval[1][0]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][1]*sqrt(pacFrac);

        for (int j=0;j<2;j++) particles[i].r[j]=actualPosition[j];
        particles[i].normR[0]=(boxMatrix[1][1]*particles[i].r[0]-boxMatrix[0][1]*particles[i].r[1])/detBoxMatrix;
        particles[i].normR[1]=-(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/detBoxMatrix;
        particles[i].phi=matrixCellPhi[cellNumber[1][0]][cellNumber[1][1]];

        columnCounter++;
        if (columnCounter>=matrixOfParticlesSize[0]*mod) {
            rowCounter++;
            columnCounter=0;
        }
    }
    if (gaps>0) return createRandomGaps(particles,boxMatrix);
    else return 1;
}

void adjustAngles (particle *particles, double boxMatrix[2][2]) {
    int tryNumber=-1,collidingPairs=0;
    printf("Angle adjusting... ");
    do {
        tryNumber++;

        collidingPairs=0;
        for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
            double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
                   normalizedRY=particles[i].normR[1]-particles[j].normR[1],
                   rx=particles[i].r[0]-particles[j].r[0],
                   ry=particles[i].r[1]-particles[j].r[1];
            rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
            ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
            double r2=rx*rx+ry*ry, dr=sqrt(r2);
            int energy;
            if (dr<maxDistance) {
                if (dr<minDistance) energy=1;
                else {
                    double gamma=atan(ry/rx),
                            aAngle=particles[j].phi-gamma,
                            bAngle=particles[i].phi-gamma;
                    if (multimerN%2!=0) {
                        if (rx>0) bAngle-=C;
                        else aAngle-=C;
                    }
                    aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                    energy=checkOverlaps(dr,aAngle,bAngle);
                }
                if (energy==1) {
                    collidingPairs++;
                    for (int k=0;k<activeN;k++) {
                        particles[k].phi+=0.00001;
                    }
                    i=activeN; j=activeN; break;
                }
            }
        }
    } while (collidingPairs>0);
    printf("Adjusted after: %d approach.\n",tryNumber);
}

void checkSinglePeriodicBoundaryConditions (particle *particle, double boxMatrix[2][2]) {
    for (int j=0;j<2;j++) {
        int change;
        do {
            change=0;
            if (particle->normR[j]<0) {
                particle->normR[j]++;
                particle->r[0]+=boxMatrix[0][j];
                particle->r[1]+=boxMatrix[1][j];
                change=1;
            } else if (particle->normR[j]>=1) {
                particle->normR[j]--;
                particle->r[0]-=boxMatrix[0][j];
                particle->r[1]-=boxMatrix[1][j];
                change=1;
            }
        } while (change==1);
    }
    //particle->phi=normalizeAngle(particle->phi);    //rotationsAnalysis(comment) 1/1
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

int getEnergy (particle *particles, int index, double boxMatrix[2][2]) {
    int energy=0;

    //matrix 4/16
    /*for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) {
        int cellIndex[2]={particles[index].cell[0]+i,particles[index].cell[1]+j};
        for (int k=0;k<2;k++) {
            if (cellIndex[k]<0) cellIndex[k]+=(int)bCSize[k];
            if (cellIndex[k]>=bCSize[k]) cellIndex[k]-=(int)bCSize[k];
        }
        for (int k=0;k<boxCellsN[cellIndex[0]][cellIndex[1]];k++) if (boxCellsIndex[cellIndex[0]][cellIndex[1]][k]!=index) {
            double normalizedRX=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].normR[0]-particles[index].normR[0],
                   normalizedRY=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].normR[1]-particles[index].normR[1],
                   rx=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].r[0]-particles[index].r[0],
                   ry=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].r[1]-particles[index].r[1];
            rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
            ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
            double r2=rx*rx+ry*ry, dr=sqrt(r2);
            if (dr<maxDistance) {
                double gamma=atan(ry/rx),
                       aAngle=particles[index].phi-gamma,
                       bAngle=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].phi-gamma;
                if (multimerN%2!=0) {
                    if (rx>0) bAngle-=C;
                    else aAngle-=C;
                }
                aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                energy=(double)checkOverlaps(dr,aAngle,bAngle);
                if (energy==1) {k=boxCellsN[cellIndex[0]][cellIndex[1]];i=2;j=2;}
            }
        }
    }*/

    for (int i=0;i<particles[index].neighCounter;i++) {
        double normalizedRX=particles[particles[index].neighbours[i]].normR[0]-particles[index].normR[0],
               normalizedRY=particles[particles[index].neighbours[i]].normR[1]-particles[index].normR[1],
               rx=particles[particles[index].neighbours[i]].r[0]-particles[index].r[0],
               ry=particles[particles[index].neighbours[i]].r[1]-particles[index].r[1];
        rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
        ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
        double r2=rx*rx+ry*ry, dr=sqrt(r2);
        if (dr<maxDistance) {
            //if (dr<minDistance) energy=1; //analyticMethodForHCM 1/8
            //else {                        //
                double gamma=atan(ry/rx),
                       aAngle=particles[index].phi-gamma,
                       bAngle=particles[particles[index].neighbours[i]].phi-gamma;
                if (multimerN%2!=0) {
                    //rozważanie, która molekuła jest 'po lewej', a która 'po prawej' (moja metoda-odwraca cząstkę po prawej; analityczna metoda-istotna kolejność)
                    //if (rx>0) bAngle-=C;    //analytiCMethodForHCM 2/8
                    //else aAngle-=C;         //
                    if (rx<0) {
                        double buffer=aAngle; aAngle=bAngle; bAngle=buffer;
                    }
                }
                //aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle); //analyticMethodForHCM 3/8
                //energy=checkOverlaps(dr,aAngle,bAngle);                       //
                energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
            //} //analyticMethodForHCM 4/8
            if (energy==1) i=particles[index].neighCounter;
        }
    }
    return energy;
}

int attemptToDisplaceAParticle (particle *particles, int index, double boxMatrix[2][2]) {
    int result=1;
    double oldR[2]={particles[index].r[0],particles[index].r[1]},
           oldNormR[2]={particles[index].normR[0],particles[index].normR[1]},
           oldPhi=particles[index].phi;
    //int oldCell[2]={particles[index].cell[0],particles[index].cell[1]}; //matrix 5/16
    for (int i=0;i<2;i++) particles[index].r[i]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaR;
    particles[index].normR[0]=(boxMatrix[1][1]*particles[index].r[0]-boxMatrix[0][1]*particles[index].r[1])/detBoxMatrix;
    particles[index].normR[1]=-(boxMatrix[1][0]*particles[index].r[0]-boxMatrix[0][0]*particles[index].r[1])/detBoxMatrix;
    particles[index].phi+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaPhi;
    while (particles[index].phi>C || particles[index].phi<C/3.0) {  //utrzymywanie kata w granicy od C/3 do C, zeby czastki sie ladnie ustawialy
        particles[index].phi=oldPhi;
        particles[index].phi+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaPhi;
    }
    //checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix); //matrix 6/16
    //for (int i=0;i<2;i++) particles[index].cell[i]=(int)(particles[index].normR[i]*bCSize[i]); //matrix 7/16
    int newEnPot=getEnergy(particles,index,boxMatrix);
    if (newEnPot==1) {
        for (int i=0;i<2;i++) {
            particles[index].r[i]=oldR[i];
            particles[index].normR[i]=oldNormR[i];
            //particles[index].cell[i]=oldCell[i]; //matrix 8/16
        }
        particles[index].phi=oldPhi;
        result=0;
    } else {
        checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix); //matrix(comment) 9/16
        /*if (particles[index].cell[0]!=oldCell[0] || particles[index].cell[1]!=oldCell[1]) { //matrix 10/16
            for (int i=0;i<boxCellsN[oldCell[0]][oldCell[1]];i++) if (boxCellsIndex[oldCell[0]][oldCell[1]][i]==index) {
                boxCellsIndex[oldCell[0]][oldCell[1]][i]=boxCellsIndex[oldCell[0]][oldCell[1]][--boxCellsN[oldCell[0]][oldCell[1]]];
                i=boxCellsN[oldCell[0]][oldCell[1]];
            }
            boxCellsIndex[particles[index].cell[0]][particles[index].cell[1]][boxCellsN[particles[index].cell[0]][particles[index].cell[1]]++]=index;
        }*/
    }
    return result;
}

int attemptToChangeVolume (particle *particles, double pressure, double boxMatrix[2][2], double *volume) {
    int result=1;
    double /*lnNewBoxMatrix[2][2], */newBoxMatrix[2][2];
    if ((*volume)/VcpPerParticle/N<pressureRealOfNotFluid) {
    //if (pressure>pressureRealOfNotFluid) {  //dozwolone zmiany ksztaltu pudla (faza stala)
        /*for (int i=0;i<2;i++) for (int j=0;j<2;j++) lnNewBoxMatrix[i][j]=log(boxMatrix[i][j]);
        lnNewBoxMatrix[0][0]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        lnNewBoxMatrix[1][1]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        if (boxMatrix[0][1]==0) lnNewBoxMatrix[0][1]=-20.0+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV; //log(0) -> -\infty; E^~0 za duze jak na start
        else lnNewBoxMatrix[0][1]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        lnNewBoxMatrix[1][0]=lnNewBoxMatrix[0][1];
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) newBoxMatrix[i][j]=exp(lnNewBoxMatrix[i][j]);*/
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        newBoxMatrix[1][1]=boxMatrix[1][1]+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        newBoxMatrix[0][1]=boxMatrix[0][1]+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        newBoxMatrix[1][0]=newBoxMatrix[0][1];
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        /*lnNewBoxMatrix[0][0]=log(boxMatrix[0][0])+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        newBoxMatrix[0][0]=exp(lnNewBoxMatrix[0][0]);*/
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;

        newBoxMatrix[1][1]=boxMatrix[1][1]/boxMatrix[0][0]*newBoxMatrix[0][0];
        newBoxMatrix[0][1]=boxMatrix[0][1]; newBoxMatrix[1][0]=boxMatrix[1][0];
    }
    double newDetBoxMatrix=newBoxMatrix[0][0]*newBoxMatrix[1][1]-newBoxMatrix[1][0]*newBoxMatrix[0][1],
           newVolume=fabs(newDetBoxMatrix);

    double bufferR[activeN][2];
    for (int i=0;i<activeN;i++) {
        bufferR[i][0]=newBoxMatrix[0][0]*particles[i].normR[0]+newBoxMatrix[0][1]*particles[i].normR[1];
        bufferR[i][1]=newBoxMatrix[1][0]*particles[i].normR[0]+newBoxMatrix[1][1]*particles[i].normR[1];
    }

    //matrix 11/16
//UWAGA: tutaj jest to zrobione na iterowaniu po komorkach, wyszlo jakies skomplikowane skanowanie tylko
//komorek z gory i z prawej, zeby nie powtarzac, jakies rozpoznawanie komorki 'srodkowej'... pojawia sie przez to duzo if'ow, etc.
//Moze lepiej zrobic analogicznie jak w przypadku listy sasiadow? Tzn. iterowac PO CZASTKACH, kazda sprawdzac z sasiadami z przyleglych
//komorek, a zeby nie bylo powtarzania, zastosowac ten sam prosty warunek: porownywac tylko, gdy index rozpatrywanej czastki jest mniejszy od porownywanej.

    /*for (int i=0;i<bCSize[0];i++) for (int j=0;j<bCSize[1];j++) {
        if (boxCellsN[i][j]>0) for (int k=0;k<5;k++) {
            int cellIndex[2];
            int zakresL=boxCellsN[i][j], startM=0;
            if (k==0) {cellIndex[0]=i; cellIndex[1]=j; zakresL--;}
            else if (k==1) {cellIndex[0]=i-1; cellIndex[1]=j+1;}
            else if (k==2) {cellIndex[0]=i; cellIndex[1]=j+1;}
            else if (k==3) {cellIndex[0]=i+1; cellIndex[1]=j+1;}
            else if (k==4) {cellIndex[0]=i+1; cellIndex[1]=j;}
            for (int l=0;l<2;l++) {

//UWAGA: tu chyba jest blad (cellIndex=-1 -> cellIndex=(int)bCSize-1, tymczasem komorka O INDEKSIE (int)bcSize TEZ moze cos zawierac, bo: .cell[i]=(int)(particles[index].normR[i]*bCSize[i])
//analogicznie, przy cellIndex=(int)bCSize+1 -> cellIndex=1 (zapomina sie o cellIndex=0). TEN SAM PROBLEM JEST przy matrix 4/16.
//Nalezy: sprawdzic dzialanie tego (wyswietlac cellIndex sprawdzanej czastki a po przecinku cellIndexy identyfikowane jako sasiednie.
//prawdopodobnie dla cellIndexu 0 i (int)bCSize blednie beda identyfikowane periodyczne sasiady.
//Poprawa kodu (na szybko myslac), nalezy dodawac/odejmowac nie (int)bCSize[l] a raczej ceil(bCSize[l]), bo (int)bCSize wskazuje na INDEKS
//najwyzszej komorki (tak zostalo to zaprojektowane), natomiast dopiero ceil(bCSize) wskazuje na index+1, czyli na LICZBE komorek (w pionie/poziomie).

                if (cellIndex[l]<0) cellIndex[l]+=(int)bCSize[l];
                if (cellIndex[l]>=bCSize[l]) cellIndex[l]-=(int)bCSize[l];
            }
            for (int l=0;l<zakresL;l++) {if (k==0) startM=l;
                for (int m=startM;m<boxCellsN[cellIndex[0]][cellIndex[1]];m++) {
                    double normalizedRX=particles[boxCellsIndex[i][j][l]].normR[0]-particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]].normR[0],
                           normalizedRY=particles[boxCellsIndex[i][j][l]].normR[1]-particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]].normR[1],
                           rx=bufferR[boxCellsIndex[i][j][l]][0]-bufferR[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]][0],
                           ry=bufferR[boxCellsIndex[i][j][l]][1]-bufferR[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]][1];
                    rx-=round(normalizedRX)*newBoxMatrix[0][0]+round(normalizedRY)*newBoxMatrix[0][1];
                    ry-=round(normalizedRX)*newBoxMatrix[1][0]+round(normalizedRY)*newBoxMatrix[1][1];
                    double r2=rx*rx+ry*ry, dr=sqrt(r2);
                    if (dr<maxDistance) {
                        double gamma=atan(ry/rx),
                               aAngle=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]].phi-gamma,
                               bAngle=particles[boxCellsIndex[i][j][l]].phi-gamma;
                        if (multimerN%2!=0) {
                            if (rx>0) bAngle-=C;
                            else aAngle-=C;
                        }
                        aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                        double energy=(double)checkOverlaps(dr,aAngle,bAngle);
                        if (energy==1) {
                            result=0;
                            l=zakresL; k=5; j=bCSize[1]; i=bCSize[0]; break;
                        }
                    }
                }
            }
        }
    }*/

    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
        if (i<particles[i].neighbours[j]) {
            double normalizedRX=particles[i].normR[0]-particles[particles[i].neighbours[j]].normR[0],
                   normalizedRY=particles[i].normR[1]-particles[particles[i].neighbours[j]].normR[1],
                   rx=bufferR[i][0]-bufferR[particles[i].neighbours[j]][0],
                   ry=bufferR[i][1]-bufferR[particles[i].neighbours[j]][1];
            rx-=round(normalizedRX)*newBoxMatrix[0][0]+round(normalizedRY)*newBoxMatrix[0][1];
            ry-=round(normalizedRX)*newBoxMatrix[1][0]+round(normalizedRY)*newBoxMatrix[1][1];
            double r2=rx*rx+ry*ry, dr=sqrt(r2);
            if (dr<maxDistance) {
                int energy;
                //if (dr<minDistance) energy=1; //analyticMethodForHCM 5/8
                //else {                        //
                    double gamma=atan(ry/rx),
                           aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                           bAngle=particles[i].phi-gamma;
                    if (multimerN%2!=0) {
                        //rozważanie, która molekuła jest 'po lewej', a która 'po prawej' (moja metoda-odwraca cząstkę po prawej; analityczna metoda-istotna kolejność)
                        //if (rx>0) bAngle-=C;    //analytiCMethodForHCM 6/8
                        //else aAngle-=C;         //
                        if (rx<0) {
                            double buffer=aAngle; aAngle=bAngle; bAngle=buffer;
                        }
                    }
                    //aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle); //analyticMethodForHCM 7/8
                    //energy=checkOverlaps(dr,aAngle,bAngle);                       //
                    energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
                //} //analyticMethodForHCM 8/8
                if (energy==1) {
                    result=0;
                    i=activeN; break;
                }
            }
        }
    }

    if (result) {
        //double arg=-(pressure*(newVolume-(*volume))-(((double)N+1.0)*log(newVolume/(*volume))+log((newBoxMatrix[0][0]+newBoxMatrix[1][1])/(boxMatrix[0][0]+boxMatrix[1][1]))));
        double arg=-(pressure*(newVolume-(*volume))+((double)N+1.0)*log(newVolume/(*volume))+log((newBoxMatrix[0][0]+newBoxMatrix[1][1])/(boxMatrix[0][0]+boxMatrix[1][1])));
        if (MTGenerate(randomStartStep)%1000000/1000000.0>exp(arg)) result=0;
        if (result) {
            *volume=newVolume;
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
            for (int i=0;i<activeN;i++) for (int j=0;j<2;j++) particles[i].r[j]=bufferR[i][j];
            detBoxMatrix=newDetBoxMatrix;
        }
    }

    return result;
}
/////////////////  } PARTICLE functions

int createIterationTable () {
    char startArguments[50]; FILE *fileStartArguments = fopen("startArguments.txt","rt");
    if (fileStartArguments==NULL) {printf("Missing file: startArguments.txt\n"); return 1;}
    while (fgets(startArguments,50,fileStartArguments)!=NULL) {
        sscanf(startArguments,"%c",startArguments); char *pEnd;
        iterationTable[fileIterateIterationsNumber][0]=strtod(startArguments,&pEnd);
        iterationTable[fileIterateIterationsNumber][1]=strtod(pEnd,&pEnd);
        iterationTable[fileIterateIterationsNumber++][2]=strtod(pEnd,NULL);
    }
    fclose(fileStartArguments);
    return 0;
}

void addAppendix (char *fileName, char *JOBID, bool jobIdOn) {
    strcpy(buffer,"2D_N-"); strncat(buffer,bufferN,20);
    strncat(buffer,"_gaps-",10); strncat(buffer,bufferGaps,20);
    strncat(buffer,"_G-",5); strncat(buffer,bufferG,5);
    strncat(buffer,"_badanie-",10); strncat(buffer,bufferFolderIndex,5);
    strncat(buffer,"_mN-",5); strncat(buffer,bufferMN,20);
    strncat(buffer,"_mS-",5); strncat(buffer,bufferMS,20);
    strncat(buffer,"_mD-",5); strncat(buffer,bufferMD,20);
    mkdir(buffer,S_IRWXU);
    strncat(buffer,"/",2);
    if (jobIdOn) {
        strncat(buffer,JOBID,50);
        strncat(buffer,"_",2);
    }
    strncat(buffer,fileName,200);
    strcpy(fileName,buffer);
}

void adjustOrientationsFile (FILE *file, char *path) {
    if (file==NULL) {
        file=fopen(path,"a"); fprintf(file,"{"); fclose(file);
    } else if (!onlyMath[0]) {
        char bufferForEraseLastChar[200],linia[110]; strcpy(bufferForEraseLastChar,path); strncat(bufferForEraseLastChar,"_BUFF",6);
        FILE *bFELC = fopen(bufferForEraseLastChar,"w");
        int poziomNawiasu=0;
        while (fgets(linia,100,file)!=NULL) {
            sscanf(linia,"%c",linia); int lastIndex=100;
            for (int i=0;i<100;i++) {
                if (linia[i]=='{') poziomNawiasu++;
                else if (linia[i]=='}') poziomNawiasu--;
                if (poziomNawiasu==0) {
                    lastIndex=i;
                    break;
                }
            }
            if (lastIndex==100) fprintf(bFELC,"%s",linia);
            else {
                for (int i=0;i<lastIndex;i++) fprintf(bFELC,"%c",linia[i]);
                fprintf(bFELC,"%c",',');
            }
        }
        fclose(bFELC); fclose(file);
        remove(path); rename(bufferForEraseLastChar,path);
    }
}

void getNextArgument (double prevArg[2], bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            for (int i=0;i<2;i++) prevArg[i]=growing?iterationTable[actIteration][i+1]:iterationTable[fileIterateIterationsNumber-1-actIteration][i+1];
            if (growing) startMinPacFrac=iterationTable[actIteration][0];
            else startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg[0]*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>=round(intervalMin[i]*10000) && round(prevArg[0]*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)+round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg[0]/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>round(intervalMin[i]*10000) && round(prevArg[0]*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)-round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)<round(minArg*10000)) growing=-1;
    }
}

double getAvErrorFromSumEps (double sum, double denominator) {
    return sqrt(sum/denominator);
}

int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[300];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,500,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20];
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++];
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:multimerN=strtol(data,NULL,10);break;
                case 4:initMode=strtol(data,NULL,10);break;
                case 5:multimerS=strtod(data,NULL);break;
                case 6:multimerD=strtod(data,NULL);break;
                case 7:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 8:growing=strtol(data,NULL,10);break;
                case 9:loadedConfiguration=strtol(data,NULL,10);break;
                case 10:loadedArg=strtod(data,NULL);break;
                case 11:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 12:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 13:loadedSetGenerator=strtol(data,NULL,10);break;
                case 14:iterationsNumber=strtol(data,NULL,10);break;
                case 15:countCollidingPairs=strtol(data,NULL,10);break;
                case 16:intervalSampling=strtol(data,NULL,10);break;
                case 17:intervalOutput=strtol(data,NULL,10);break;
                case 18:saveConfigurations=strtol(data,NULL,10);break;
                case 19:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 20:ODFLength=strtol(data,NULL,10);break;
                case 21:OCFMode=strtol(data,NULL,10);break;
                case 22:neighRadiusMod=strtod(data,NULL);break;
                case 23:intervalOrientations=strtol(data,NULL,10);break;
                case 24:skipFirstIteration=strtol(data,NULL,10);break;
                case 25:useSpecificDirectory=strtol(data,NULL,10);break;
                case 26:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 27:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 28:intervalResults=strtol(data,NULL,10);break;
                case 29:maxDeltaR=strtod(data,NULL);break;
                case 30:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 31:desiredAcceptanceRatioV=strtod(data,NULL);break;
                case 32:useFileToIterate=strtol(data,NULL,10);break;
                case 33:startMinPacFrac=strtod(data,NULL);break;
                case 34:startMaxPacFrac=strtod(data,NULL);break;
                case 35:minArg=strtod(data,NULL);break;
                case 36:maxArg=strtod(data,NULL);break;
                case 37:multiplyArgument=strtol(data,NULL,10);break;
                case 38:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-39)%3) {
                        case 0: intervalMin[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 1: intervalMax[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 2: intervalDelta[intervalLicznik++/3]=strtod(data,NULL);break;
                    }
                    break;
            }
            dataIndex++;
            for (int i=0;i<20;i++) data[i]=' ';
        }
        fclose(fileConfig);
    } while (testValue!=12345);

    //zlecanie parametrow z poziomu wiersza polecen:
    char JOBID[50]="j-"; int pointNumber=0;
    if (argumentsNumber==1) {
        strncat(JOBID,"none",50);
        if (useFileToIterate) if(createIterationTable()) return 0;
    } else {
        int correctNumberOfArguments=1;
        switch (strtol(arguments[1],NULL,10)) {
            case 0: //ustaw JOBID
                if (argumentsNumber==3) {
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (useFileToIterate) if(createIterationTable()) return 0;
                } else correctNumberOfArguments=0; break;
            case 1: //ustaw JOBID, singleRun dla parametrow zadanych bezposrednio
                if (argumentsNumber==13) {
                    useFileToIterate=0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (growing) {
                        startMinPacFrac=strtod(arguments[3],NULL); minArg=strtod(arguments[4],NULL);
                    } else {
                        startMaxPacFrac=strtod(arguments[3],NULL); maxArg=strtod(arguments[4],NULL);
                    }
                    N=strtol(arguments[5],NULL,10);
                    gaps=strtol(arguments[6],NULL,10);
                    multimerS=strtod(arguments[7],NULL);
                    multimerD=strtod(arguments[8],NULL);
                    growing=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    multimerN=strtol(arguments[12],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    N=strtol(arguments[3],NULL,10);
                    gaps=strtol(arguments[4],NULL,10);
                    multimerS=strtod(arguments[5],NULL);
                    multimerD=strtod(arguments[6],NULL);
                    growing=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    pointNumber=strtol(arguments[11],NULL,10);
                    generatorStartPoint=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1;
                    loadedArg=strtod(arguments[9],NULL);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1; loadType=1;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==13) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    onlyMath[0]=1;
                    onlyMath[1]=strtol(arguments[3],NULL,10);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=0;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    multimerN=strtol(arguments[12],NULL,10);
                    skipFirstIteration=0;
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, multimerN\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, generatorStartPoint, multimerN\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration, multimerN\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, multimerN\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, multimerN\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        if (growing) {
            startMinPacFrac=iterationTable[0][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[0][i+1];
        } else {
            startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[fileIterateIterationsNumber-1][i+1];
        }
    } else {
        startArg[0]=growing?minArg:maxArg;
        startArg[1]=0;
    }
    pressureRealOfNotFluid/=(multimerS*multimerS);
    deltaR=maxDeltaR*multimerS; deltaV*=multimerS;
    for (int i=0;i<pointNumber;i++) getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg[0];
    activeN=N-gaps;

    if (N%56!=0 && N%780!=0 && floor(sqrt(N))!=sqrt(N)) {
        printf("ERROR: Not supported N: %d.\n",N);
        return 0;
    }

    //stale wynikajace z zadanych parametrow multimerow
    L=multimerS/multimerD;
    C=pi/(double)multimerN; deltaPhi=deltaR*2.0*sin(C);
    ROkreguOpisanego=multimerS/(2.0*sin(C));
    absoluteMinimum=atan(L/(2.0*L/tan(C)+sqrt(4.0-L*L)));
    absoluteMinimum2=C-absoluteMinimum;
    if (multimerN%2==0) minDistance=getMinimalDistanceAnalyticalMethodForEvenHCM(absoluteMinimum,absoluteMinimum);
    else minDistance=getMinimalDistanceAnalyticalMethodForOddHCM(absoluteMinimum,absoluteMinimum-C);
    maxDistance=ROkreguOpisanego*2+multimerD;
    neighRadius=neighRadiusMod*maxDistance; neighRadius2=neighRadius*neighRadius; neighSafeDistance=neighRadius-maxDistance;

    //nazwy folderow na podstawie parametrow programu
    sprintf(bufferG,"%d",growing); sprintf(bufferN,"%d",N); sprintf(bufferGaps,"%d",gaps);
    sprintf(bufferMN,"%d",multimerN); sprintf(bufferMS,"%.2f",multimerS); sprintf(bufferMD,"%.6f",multimerD);

    int folderIndex=useSpecificDirectory, checkNext;
    char bufferCheckFolderExisting[200];
    FILE *checkFolderExisting;
    if (!folderIndex) do {
        sprintf(bufferFolderIndex,"%d",++folderIndex);
        strcpy(bufferCheckFolderExisting,resultsFileName);
        addAppendix(bufferCheckFolderExisting,JOBID,false);
        checkFolderExisting = fopen(bufferCheckFolderExisting,"rt");
        if (checkFolderExisting!=NULL) {
            fclose(checkFolderExisting);
            checkNext=1;
        } else checkNext=0;
    } while (checkNext);
    sprintf(bufferFolderIndex,"%d",folderIndex);
    addAppendix(resultsFileName,JOBID,false);
    addAppendix(excelResultsFileName,JOBID,false);
    addAppendix(configurationsFileName,JOBID,true);
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.3f",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(orientationsFileName,JOBID,true);
    if (OCFMode) addAppendix(orientatCorrelFunFileName,JOBID,true);
    addAppendix(orientationsResultsFileName,JOBID,true);
    addAppendix(configurationsListFileName,JOBID,false);

    particle particles[N];
    double args[10];

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileOrientations, *fileOrientatCorrelFun, *fileConfigurationsList, *fileAllResults, *fileAllOrientations, *fileOrientationsResults, *fileAllOrientationsResults;
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        if (saveConfigurations) fprintf(fileResults,"Cycles\tPressureReduced\tVolume\tBoxMatrix[0][0]\tBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tRho\tV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavB\tdAvB\tavMy\tdAvMy\tavE\tdAvE\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileResults,"Cycles\tPressureReduced\tVolume\tBoxMatrix[0][0]\tBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tRho\tV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavB\tdAvB\tavMy\tdAvMy\tavE\tdAvE\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        if (saveConfigurations) fprintf(fileExcelResults,"PressureReduced\tV/V_cp\tavNu\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tavB\tavMy\tavE\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileExcelResults,"PressureReduced\tV/V_cp\tavNu\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tavB\tavMy\tavE\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg[2]={startArg[0],startArg[1]}, oldBoxMatrix[2][2];
    while (growing>=0) {
        double pressureReduced=arg[0], pressureReal=pressureReduced/multimerS/multimerS, //pressureReduced=\tau*\sigma^2/(kT), kT=1
               boxMatrix[2][2],matrixOfParticlesSize[2],unitCellAtCP[2],
               matrixCellXY[6][6][2],matrixCellPhi[6][6];
        int n[2]; //n[X/Y], matrixCell[n[XMax]][n[YMax]][x/y], zatem: n[X/Y](max)=6
        switch (multimerN) {
            /*case 6: {//dla heksamerow o dowolnym d/\sigma
                unitCellAtCP[0]=minDistance; unitCellAtCP[1]=sqrt(3)*minDistance;
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=absoluteMinimum2;
                matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=absoluteMinimum2;
            } break;*/
            case 5: {//dla pentamerow o d/\sigma=1
                unitCellAtCP[0]=2.4048671732*multimerS; unitCellAtCP[1]=4.2360679772*multimerS;
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                matrixCellXY[0][1][0]=1.0131106571*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
            } break;
            case 7: switch (initMode) {//dla heptamerow o d/\sigma=1, struktury jak w WojTreKow2003PRE
                case 0: {//struktura A
                    unitCellAtCP[0]=3.0566685376*multimerS; unitCellAtCP[1]=5.5150210832*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=1.9144263193*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//struktura B
                    unitCellAtCP[0]=3.0566685376*multimerS; unitCellAtCP[1]=5.5382990167*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=1.3656999121*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 2: {//struktura C
                    unitCellAtCP[0]=6.0309063912*multimerS; unitCellAtCP[1]=5.5371391576*multimerS;
                    n[0]=2; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.1644029030;
                    matrixCellXY[1][0][0]=unitCellAtCP[0]/2.0; matrixCellXY[1][0][1]=0.1230590461*multimerS; matrixCellPhi[1][0]=0.1644029030;
                    matrixCellXY[0][1][0]=1.2832327269*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.2843960475;
                    matrixCellXY[1][1][0]=4.2986859225*multimerS; matrixCellXY[1][1][1]=2.8916286250*multimerS; matrixCellPhi[1][1]=0.2843960475;
                } break;
            } break;
            /*case 95: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=31.1028*multimerS; unitCellAtCP[1]=54.084*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.5404*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie (przechodzi w DISORDERED nawet przy p*=4000)
                    unitCellAtCP[0]=31.1028*multimerS; unitCellAtCP[1]=54.1536*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.5423*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 97: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=31.7394*multimerS; unitCellAtCP[1]=55.1314*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.8679*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie
                    unitCellAtCP[0]=31.7394*multimerS; unitCellAtCP[1]=55.2852*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.8679*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 98: switch (initMode) {
                case 0: {//naprzemienne katy gamma (jednakowe gamma w rzedzie)
                    unitCellAtCP[0]=32.0695*multimerS; unitCellAtCP[1]=55.7368*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=0.0155995;
                    matrixCellXY[0][1][0]=16.0347*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.0485146;
                } break;
                case 1: {//jednakowe gamma wszedzie
                    unitCellAtCP[0]=32.0695*multimerS; unitCellAtCP[1]=55.7804*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=0.01559950076973243;
                    matrixCellXY[0][1][0]=16*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.01559950076973243;
                } break;
            } break;
            case 100: switch (initMode) {
                case 0: {//naprzemienne katy gamma (jednakowe gamma w rzedzie)
                    unitCellAtCP[0]=32.6904*multimerS; unitCellAtCP[1]=56.7792*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0161203;
                    matrixCellXY[0][1][0]=16.3452*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.0161203;
                } break;
                case 1: {//jednakowe gamma wszedzie
                    unitCellAtCP[0]=32.6904*multimerS; unitCellAtCP[1]=56.9516*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0161203;
                    matrixCellXY[0][1][0]=16.3452*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.0161203;
                } break;
            } break;
            case 101: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=33.0128*multimerS; unitCellAtCP[1]=57.4142*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.5064*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie
                    unitCellAtCP[0]=33.0128*multimerS; unitCellAtCP[1]=57.435*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.5064*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 103: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=33.6495*multimerS; unitCellAtCP[1]=58.4394*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.8248*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie
                    unitCellAtCP[0]=33.6495*multimerS; unitCellAtCP[1]=58.5904*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.8248*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 104: switch (initMode) {
                case 0: {//naprzemienne katy gamma (jednakowe gamma w rzedzie)
                    unitCellAtCP[0]=33.963916981946*multimerS; unitCellAtCP[1]=59.062029530994*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0154856;
                    matrixCellXY[0][1][0]=16.982008490998*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.0154856;
                } break;
                case 1: {//jednakowe gamma wszedzie
                    unitCellAtCP[0]=33.963916981946*multimerS; unitCellAtCP[1]=59.104*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0154856;
                    matrixCellXY[0][1][0]=16.982008490998*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.0154856;
                } break;
            } break;*/
            default: {
                unitCellAtCP[0]=maxDistance; unitCellAtCP[1]=maxDistance*sqrt(3);
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
            } break;
        } VcpPerParticle=unitCellAtCP[0]*unitCellAtCP[1]/(double)n[0]/(double)n[1];
        if (N%56==0) {matrixOfParticlesSize[0]=7; matrixOfParticlesSize[1]=8;}
        else if (N%780==0) {matrixOfParticlesSize[0]=26; matrixOfParticlesSize[1]=30;}
        else if (floor(sqrt(N))==sqrt(N)) {matrixOfParticlesSize[0]=matrixOfParticlesSize[1]=sqrt(N);}
        double NLinearMod = sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]);
        if (startArg[0]==arg[0]) {
            for (int i=0;i<2;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]/(double)n[i]*NLinearMod*(growing?sqrt(startMinPacFrac):sqrt(startMaxPacFrac));
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0;
        } else for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
        detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
        double volume=fabs(detBoxMatrix), rho=N/volume, pacFrac=1.0/VcpPerParticle/rho;

        if (!onlyMath[0]) {
            if (arg[0]==startArg[0] && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (StartDen: %.7E, startPacFrac: %.7E), mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,startArg[0],rho,pacFrac,multimerN,multimerS,multimerD);
                if (!initPositions(particles,boxMatrix,matrixOfParticlesSize,n,matrixCellXY,matrixCellPhi,pacFrac)) return 0;
                adjustAngles(particles,boxMatrix);
                updateNeighbourList(particles,boxMatrix); //matrix(comment) 12/16
            } else if (loadedConfiguration) {
                char configurations[4096];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                for (int i=0;i<3;i++) fgets(configurations,4096,fileCTL);
                int character,dataType=0,pIndex=-1; char data[50]=""; int actIndex=0;
                while (dataType<11) {
                    character=fgetc(fileCTL); //character is in int, but it can work as char
                    if (dataType<10) { //stage #1 configuration parameters
                        if (character==' ') {
                            data[actIndex++]=' '; //end of data without clearing the entire array
                            args[dataType++]=strtod(data,NULL);
                            if (dataType==10) {
                                boxMatrix[0][0]=args[5]; boxMatrix[1][1]=args[6]; boxMatrix[1][0]=args[7]; boxMatrix[0][1]=args[7];
                                deltaR=args[8]; deltaPhi=deltaR*2.0*sin(C); deltaV=args[9];
                                arg[0]=args[3]; pressureReduced=arg[0];

                                detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
                                volume=fabs(detBoxMatrix); rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            }
                            actIndex=0;
                            continue;
                        }
                        data[actIndex++]=character;
                    } else { //stage #2 configuration parameters (coordinates of particles)
                        if (character=='m') {pIndex++; continue;}
                        if (pIndex>=0) {
                            for (int i=0;i<3;i++) {
                                character=fgetc(fileCTL);
                                while (character!=',') {
                                    data[actIndex++]=character;
                                    character=fgetc(fileCTL);
                                } data[actIndex++]=' '; actIndex=0;
                                if (i<2) particles[pIndex].r[i]=strtod(data,NULL);
                                else particles[pIndex].phi=strtod(data,NULL);
                            }
                            particles[pIndex].normR[0]=(boxMatrix[1][1]*particles[pIndex].r[0]-boxMatrix[0][1]*particles[pIndex].r[1])/detBoxMatrix;
                            particles[pIndex].normR[1]=-(boxMatrix[1][0]*particles[pIndex].r[0]-boxMatrix[0][0]*particles[pIndex].r[1])/detBoxMatrix;
                            while (character!=']') character=fgetc(fileCTL); fgetc(fileCTL); //next to read: 'm'
                            if (pIndex>=activeN-1) dataType++;
                        }
                    }
                }
                fclose(fileCTL);
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (startDen: %.7E, startPacFrac: %.7E), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4E, DeltaV: %.4E\n",N,gaps,growing,args[3],args[2],pacFrac,args[0],args[1],(long)args[4],args[8],args[9]);
                //for (int i=0;i<2;i++) for (int j=0;j<2;j++) printf("boxMatrix[%d][%d]=%.17E\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<activeN;i++) printf("%d: %.17E,  %.17E,  %.17E\n",i,particles[i].r[0],particles[i].r[1],particles[i].phi);return 0;
                updateNeighbourList(particles,boxMatrix); //matrix(comment) 13/16
            }
        }

        if (skipFirstIteration) {
            printf("Skipping first iteration...\n");
            skipFirstIteration=0;
        } else {
            if (!onlyMath[0]) {
                if (loadedConfiguration) {
                    if (loadedSetStartGenerator) {
                        printf("Setting start position of p-random number generator to position from file...\n");
                        InitMT((unsigned int)args[0]);
                        randomStartStep[0]=args[0];
                    } else {
                        randomStartStep[0]=InitRandomMT();
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                    }
                    randomStartStep[1]=0;
                    if (loadedSetGenerator) {
                        printf("Setting p-random number generator to last position from file...\n");
                        for (double i=0;i<args[1];i++) MTGenerate(randomStartStep);
                    } else printf("Setting p-random number generator to last position from file - DISABLED\n");
                } else {
                    if (generatorStartPoint==0) {
                        generatorStartPoint=time(0);
                        printf("Setting start position of p-random number generator to actual CPU time...\n");
                    } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                    InitMT((unsigned int)generatorStartPoint);
                    randomStartStep[0]=generatorStartPoint;
                    randomStartStep[1]=0;
                }
                printf("Start of equilibration at reduced pressure: %.7E (startDen: %.7E, startPacFrac: %.7E)... (%ld cycles)\n",pressureReduced,rho,pacFrac,cyclesOfEquilibration);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.7E, mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,pressureReduced,multimerN,multimerS,multimerD);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                fullCycle=activeN+volumeMoveChance,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberV=0, displacedNumberV=0,
                cyclesOfMeasurementBuffer=arg[1]==0?cyclesOfMeasurement:0;
            double nStep=fullCycle*((double)cyclesOfEquilibration+(double)cyclesOfMeasurementBuffer+10.0),
                   possibleDistance=0; //matrix(comment) 14/16
            int volumeMove=0, cycleCounter=0, indexScanned=(matrixOfParticlesSize[0]*round(matrixOfParticlesSize[1]*NLinearMod/2.0)-round(matrixOfParticlesSize[0]/2.0))*NLinearMod;
            //assignParticlesToCells(particles,boxMatrix); //matrix 15/16

            char allResultsFileName[200],bufferConfigurations[200],bufferSavedConfigurations[200],bufferOrientations[200],bufferOrientatCorrelFun[200],allOrientationsFileName[200],bufferOrientationsResults[200],allOrientationsResultsFileName[200],bufferPressure[100];
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferOrientations,orientationsFileName); strcpy(allOrientationsFileName,orientationsFileName);
            if (OCFMode) strcpy(bufferOrientatCorrelFun,orientatCorrelFunFileName);
            strcpy(bufferOrientationsResults,orientationsResultsFileName); strcpy(allOrientationsResultsFileName,orientationsResultsFileName);
            sprintf(bufferPressure,"%.3f",pressureReduced);
            strncat(allResultsFileName,"_arg-",6); strncat(allResultsFileName,bufferPressure,100); strncat(allResultsFileName,"_Results.txt",13);
            strncat(bufferConfigurations,"_arg-",6); strncat(bufferConfigurations,bufferPressure,100); strcpy(bufferSavedConfigurations,bufferConfigurations); strncat(bufferConfigurations,".txt",5); strncat(bufferSavedConfigurations,"_transient.txt",15);
            strncat(bufferOrientations,"_arg-",6); strncat(bufferOrientations,bufferPressure,100); strncat(bufferOrientations,".txt",5);
            strncat(allOrientationsFileName,"_arg-",6); strncat(allOrientationsFileName,bufferPressure,100); strncat(allOrientationsFileName,"_allOnt.txt",12);
            if (OCFMode) {
                strncat(bufferOrientatCorrelFun,"_arg-",6); strncat(bufferOrientatCorrelFun,bufferPressure,100); strncat(bufferOrientatCorrelFun,".txt",5);
            }
            strncat(bufferOrientationsResults,"_arg-",6); strncat(bufferOrientationsResults,bufferPressure,100); strncat(bufferOrientationsResults,".txt",5);
            strncat(allOrientationsResultsFileName,"_arg-",6); strncat(allOrientationsResultsFileName,bufferPressure,100); strncat(allOrientationsResultsFileName,"_allOnt.txt",12);

            fileAllResults = fopen(allResultsFileName,"a");
            adjustOrientationsFile(fileOrientations=fopen(bufferOrientations,"rt"),bufferOrientations);
            fileOrientations = fopen(bufferOrientations,"a");
            fileAllOrientations = fopen(allOrientationsFileName,"a");
            if (saveConfigurations) fileSavedConfigurations = fopen(bufferSavedConfigurations,"a");
            if (OCFMode) {
                adjustOrientationsFile(fileOrientatCorrelFun=fopen(bufferOrientatCorrelFun,"rt"),bufferOrientatCorrelFun);
                fileOrientatCorrelFun = fopen(bufferOrientatCorrelFun,"a");
            }
            if (onlyMath[0]) nStep=0;

            timeStart=time(0);
            for (double iStep=1;iStep<nStep;iStep++) {
                int randIndex;
                if (volumeMove) {
                    randIndex = (int)(MTGenerate(randomStartStep)%1000000/1000000.0*activeN);
                    volumeMove=0;
                } else randIndex = (int)(MTGenerate(randomStartStep)%1000000/1000000.0*fullCycle);
                if (randIndex<activeN) {
                    attemptedNumberR++;
                    if (attemptToDisplaceAParticle(particles,randIndex,boxMatrix))
                        displacedNumberR++;
                } else {
                    volumeMove=1;
                    attemptedNumberV++;
                    if (attemptToChangeVolume(particles,pressureReal,boxMatrix,&volume))
                        displacedNumberV++;
                }

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++;

                    if (cycle%intervalSampling==0) {
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR,
                               acceptanceRatioV = displacedNumberV/(double)attemptedNumberV;
                        /*possibleDistance+=sqrt(0.5*deltaR*deltaR)*((double)intervalSampling)*acceptanceRatioR*5.0;  //ostatnie *5.0 - dla bezpieczenstwa; sqrt(0.5*deltaR*0.5*deltaR+0.5*deltaR*0.5*deltaR)=sqrt(0.5*deltaR*deltaR)
                        if (possibleDistance>=neighSafeDistance) {
                            updateNeighbourList(particles,boxMatrix);
                            possibleDistance=0;
                        }*/  //matrix(comment) 16/16

                        /////wypisywanie danych czesciej niz normalnie i PRZED zrownowagowaniem
                        /*if (cycle%50==0) {
                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
                                           normalizedRY=particles[i].normR[1]-particles[j].normR[1],
                                           rx=particles[i].r[0]-particles[j].r[0],
                                           ry=particles[i].r[1]-particles[j].r[1];
                                    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                                    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                                    double r2=rx*rx+ry*ry, dr=sqrt(r2);
                                    int energy;
                                    if (dr<maxDistance) {
                                        if (dr<minDistance) energy=1;
                                        else {
                                            double gamma=atan(ry/rx),
                                                    aAngle=particles[j].phi-gamma,
                                                    bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (rx>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                                            energy=checkOverlaps(dr,aAngle,bAngle);
                                            //energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
                                            if (energy==1) printf("colliding: distance- %.12E, alpha: %.12E, beta: %.12E, analytic: %.12E, i: %d, j: %d\n",dr,aAngle,bAngle,getMinimalDistanceAnalyticalMethodForEvenHCM(normalizeAngle(aAngle+C),normalizeAngle(bAngle+C)),i,j);
                                        }
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+arg5),collidingPairs);
                            else printf("Cycle: %ld\n",(cycle+arg5));
                            printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                            printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                            printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                        }*/
                        /////

                        if (cycle>cyclesOfEquilibration) {
                            if (timeEquilibration==0) timeEquilibration=time(0);

                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
                                           normalizedRY=particles[i].normR[1]-particles[j].normR[1],
                                           rx=particles[i].r[0]-particles[j].r[0],
                                           ry=particles[i].r[1]-particles[j].r[1];
                                    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                                    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                                    double r2=rx*rx+ry*ry, dr=sqrt(r2);
                                    int energy;
                                    if (dr<maxDistance) {
                                        if (dr<minDistance) energy=1;  //analyticCheck 1/3
                                        else {                         //
                                            double gamma=atan(ry/rx),
                                                    aAngle=particles[j].phi-gamma,
                                                    bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (rx>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);  //analyticCheck 2/3
                                            energy=checkOverlaps(dr,aAngle,bAngle);                        //
                                            //energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
                                        }  //analyticCheck 3/3
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0)
                                fprintf(fileAllResults,"%.17E\t%.17E\t%.17E\t\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);

                            if (cycle%intervalOrientations==0) {
                                /////skan po konkretnej cząstce (np. dla 224: 14*(16/2)-(14/2)=105, etc.) - w srodku by nie skakala na granicy pudla periodycznego; bezposrednie uzycie w Mathematice (format tablicy)
                                if (cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer)
                                    fprintf(fileOrientations,"{%.12E,%.12E,%.12E}}",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                else fprintf(fileOrientations,"{%.12E,%.12E,%.12E},",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                /////skan po wszystkich cząstkach
                                for (int i=0;i<activeN-1;i++) fprintf(fileAllOrientations,"%.12E,",particles[i].phi);
                                fprintf(fileAllOrientations,"%.12E\n",particles[activeN-1].phi);
                                /////OCF
                                if (OCFMode) {
                                    char bufferText[4096]="", bufferAngle[20];
                                    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                        if (i<particles[i].neighbours[j]) {
                                            strcpy(bufferText,"");
                                            double normalizedRX=particles[i].normR[0]-particles[particles[i].neighbours[j]].normR[0],
                                                   normalizedRY=particles[i].normR[1]-particles[particles[i].neighbours[j]].normR[1],
                                                   rx=particles[i].r[0]-particles[particles[i].neighbours[j]].r[0],
                                                   ry=particles[i].r[1]-particles[particles[i].neighbours[j]].r[1];
                                            rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                                            ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                                            double gamma=atan(ry/rx),
                                                   aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                                                   bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (rx>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle+C); bAngle=normalizeAngle(bAngle+C);

                                            strncat(bufferText,"{",2); sprintf(bufferAngle,"%.12E",aAngle); strncat(bufferText,bufferAngle,20);
                                            strncat(bufferText,",",2); sprintf(bufferAngle,"%.12E",bAngle); strncat(bufferText,bufferAngle,20);
                                            if (i>=activeN-2 && cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer) strncat(bufferText,"}}",3); // -2 because it's last particle which has neighbour UNtested (the last one)
                                            else strncat(bufferText,"},",3);
                                            fprintf(fileOrientatCorrelFun,"%s",bufferText);
                                        }
                                    }
                                }
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.12E\t%.12E\t%.12E\t{",(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++)
                                    fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E],",particles[i].r[0],particles[i].r[1],particles[i].phi);
                                fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E]}",particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].phi);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+(long)args[4]),collidingPairs);
                                else printf("Cycle: %ld\n",(cycle+(long)args[4]));
                                printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                                printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                                for (int i=0;i<activeN;i++)
                                    fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                fclose(fileConfigurations);
                            }
                        }// else {//dostosowywanie delt w trakcie pomiaru
                            if (acceptanceRatioR>desiredAcceptanceRatioR) deltaR*=1.05; else deltaR*=0.95;
                            if (deltaR>maxDeltaR) deltaR=maxDeltaR;
                            deltaPhi=deltaR*2.0*sin(C);
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; else deltaV*=0.95;
                        //}
                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberV=0; displacedNumberV=0;
                    }
                    if (cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer) iStep=nStep;
                }
            }
            fclose(fileAllResults);
            fclose(fileOrientations); fclose(fileAllOrientations);
            if (saveConfigurations) fclose(fileSavedConfigurations); if (OCFMode) fclose(fileOrientatCorrelFun);
            if (timeEquilibration==0) timeEquilibration=time(0);
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            printf("Start of calculation of results...\n");

            //obliczenie liczby linii danych (potrzebne do podziału na zespoły i obliczenia średnich błędów)
            printf("Calculation of data lines... "); fflush(stdout);
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[4096]; double dataLicznik=0;
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            while(fgets(linia,300,fileAllResults)!=NULL) dataLicznik++;
            fclose(fileAllResults);
            printf("done (found %ld data lines)\n",(long)dataLicznik);

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... "); fflush(stdout);
            double avVolumeSet[10], avBoxMatrixSet[10][3], avRhoSet[10], avPacFracSet[10]; //wyniki dzielone są na 10 zespołów (obliczane są nieskorelowane wzajemnie średnie "lokalne", do obliczenia błędu średniej "globalnej")
            for (int i=0;i<10;i++) {
                avVolumeSet[i]=0; for (int j=0;j<3;j++) avBoxMatrixSet[i][j]=0;
                avRhoSet[i]=0; avPacFracSet[i]=0;
            }
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            double lineCounter=0;
            /*while(fgets(linia,300,fileAllResults)!=NULL) {//OLD-VERSION - large files including: cycles,volume,3xBoxMatrix,rho,pacFrac  1/2
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++; if (actIndex>=300) continue;
                int dataIndex=0; double dataD[6]; while (dataIndex<6) {
                    char data[50]="";
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    dataD[dataIndex++]=strtod(data,NULL);
                }
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    avVolumeSet[setIndex]+=dataD[0]; avBoxMatrixSet[setIndex][0]+=dataD[1];
                    avBoxMatrixSet[setIndex][1]+=dataD[2]; avBoxMatrixSet[setIndex][2]+=dataD[3];
                    avRhoSet[setIndex]+=dataD[4]; avPacFracSet[setIndex]+=dataD[5];
                }
                lineCounter++;
            }
            fclose(fileAllResults);
            double avVolume=0, avBoxMatrix[3]={0,0,0}, avRho=0, avPacFrac=0;
            for (int i=0;i<10;i++) {
                avVolumeSet[i]/=dataLicznik*0.1; avVolume+=avVolumeSet[i];
                avRhoSet[i]/=dataLicznik*0.1; avRho+=avRhoSet[i];
                avPacFracSet[i]/=dataLicznik*0.1; avPacFrac+=avPacFracSet[i];
                for (int j=0;j<3;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
            }*/
            while(fgets(linia,300,fileAllResults)!=NULL) {//NEW-VERSION - small files including only: 3xBoxMatrix  1/2
                sscanf(linia,"%c",linia);
                int actIndex=0;
                int dataIndex=0; double dataD[3]; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    dataD[dataIndex++]=strtod(data,NULL);
                }
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    for (int i=0;i<3;i++) avBoxMatrixSet[setIndex][i]+=dataD[i];
                }
                lineCounter++;
            }
            fclose(fileAllResults);
            double avVolume=0, avBoxMatrix[3]={0,0,0}, avRho=0, avPacFrac=0;
            for (int i=0;i<10;i++) {
                for (int j=0;j<3;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
                avVolumeSet[i]=fabs(avBoxMatrixSet[i][0]*avBoxMatrixSet[i][1]-avBoxMatrixSet[i][2]*avBoxMatrixSet[i][2]); avVolume+=avVolumeSet[i];
                avRhoSet[i]=N/avVolumeSet[i]; avRho+=avRhoSet[i];
                avPacFracSet[i]=1.0/VcpPerParticle/avRhoSet[i]; avPacFrac+=avPacFracSet[i];
            }
            avVolume*=0.1; avRho*=0.1; avPacFrac*=0.1; for (int i=0;i<3;i++) avBoxMatrix[i]*=0.1;
            //obliczenie bledow mierzonych wielkosci
            double dAvVolume=0, dAvBoxMatrix[3]={0,0,0}, dAvRho=0, dAvPacFrac=0;
            for (int i=0;i<10;i++) {double epsilon=avVolume-avVolumeSet[i]; dAvVolume+=epsilon*epsilon;} dAvVolume=getAvErrorFromSumEps(dAvVolume,90.0); //10*9 (n(n-1))
            for (int j=0;j<3;j++) {for (int i=0;i<10;i++) {double epsilon=avBoxMatrix[j]-avBoxMatrixSet[i][j]; dAvBoxMatrix[j]+=epsilon*epsilon;} dAvBoxMatrix[j]=getAvErrorFromSumEps(dAvBoxMatrix[j],90.0);}
            for (int i=0;i<10;i++) {double epsilon=avRho-avRhoSet[i]; dAvRho+=epsilon*epsilon;} dAvRho=getAvErrorFromSumEps(dAvRho,90.0);
            for (int i=0;i<10;i++) {double epsilon=avPacFrac-avPacFracSet[i]; dAvPacFrac+=epsilon*epsilon;} dAvPacFrac=getAvErrorFromSumEps(dAvPacFrac,90.0);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... "); fflush(stdout);
            double e1111Set[10], e1122Set[10], e1212Set[10], e2222Set[10], e1112Set[10], e1222Set[10],
                   HxyHyx,HxxHyy,HxxHxy,Hxx2,Hyy2,mod0,mod1;
            for (int i=0;i<10;i++) {e1111Set[i]=0; e1122Set[i]=0; e1212Set[i]=0; e2222Set[i]=0; e1112Set[i]=0; e1222Set[i]=0;}
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            lineCounter=0; int setIndex, oldSetIndex=-1;
            while(fgets(linia,300,fileAllResults)!=NULL) {
                setIndex=(int)(lineCounter/dataLicznik*10.0);
                if (setIndex!=oldSetIndex) {
                    HxyHyx=avBoxMatrixSet[setIndex][2]*avBoxMatrixSet[setIndex][2];
                    HxxHyy=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][1];
                    HxxHxy=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][2];
                    Hxx2=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][0];
                    Hyy2=avBoxMatrixSet[setIndex][1]*avBoxMatrixSet[setIndex][1];
                    mod0=HxyHyx-HxxHyy;
                    mod1=1.0/(2.0*mod0*mod0);
                    oldSetIndex=setIndex;
                }
                sscanf(linia,"%c",linia);
                int actIndex=0;
                /*while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++;  //OLD-VERSION - large files including: cycles,volume,3xBoxMatrix,rho,pacFrac  2/2
                if (actIndex>=300) continue;
                double h11,h22,h12;
                int dataIndex=0; while (dataIndex<4) {
                    char data[50]="";
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    switch (dataIndex++) {
                        case 1: h11=strtod(data,NULL);break;
                        case 2: h22=strtod(data,NULL);break;
                        case 3: h12=strtod(data,NULL);break;
                    }
                }*/
                double h11,h22,h12;  //NEW-VERSION - small files including only: 3xBoxMatrix  2/2
                int dataIndex=0; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    switch (dataIndex++) {
                        case 0: h11=strtod(data,NULL);break;
                        case 1: h22=strtod(data,NULL);break;
                        case 2: h12=strtod(data,NULL);break;
                    }
                }
                if (dataIndex<10) {
                    double hxyhyx=h12*h12,
                           hxxPhyy=h11+h22,
                           hxx2=h11*h11,
                           hyy2=h22*h22,

                           e11=mod1*(HxyHyx*(hxyhyx-HxyHyx+hyy2)-2.0*(h11*avBoxMatrixSet[setIndex][2]*h12-avBoxMatrixSet[setIndex][0]*HxyHyx+avBoxMatrixSet[setIndex][2]*h12*h22)*avBoxMatrixSet[setIndex][1]+(hxx2-Hxx2+hxyhyx)*Hyy2),
                           e22=mod1*(HxyHyx*(hxx2+hxyhyx-HxyHyx)-2.0*(HxxHxy*h12*hxxPhyy-HxxHyy*HxyHyx)+Hxx2*(hxyhyx+hyy2-Hyy2)),
                           e12=mod1*(-HxxHxy*(hxyhyx+hyy2)+HxxHyy*h12*hxxPhyy+avBoxMatrixSet[setIndex][2]*(h12*avBoxMatrixSet[setIndex][2]*hxxPhyy-(hxx2+hxyhyx)*avBoxMatrixSet[setIndex][1]));

                    e1111Set[setIndex]+=e11*e11;
                    e1122Set[setIndex]+=e11*e22;
                    e1212Set[setIndex]+=e12*e12;
                    e2222Set[setIndex]+=e22*e22;
                    e1112Set[setIndex]+=e11*e12;
                    e1222Set[setIndex]+=e12*e22;
                }
                lineCounter++;
            }
            fclose(fileAllResults);
            double e1111=0, e1122=0, e1212=0, e2222=0, e1112=0, e1222=0;
            for (int i=0;i<10;i++) {
                e1111Set[i]/=dataLicznik*0.1; e1111+=e1111Set[i];
                e1122Set[i]/=dataLicznik*0.1; e1122+=e1122Set[i];
                e1212Set[i]/=dataLicznik*0.1; e1212+=e1212Set[i];
                e2222Set[i]/=dataLicznik*0.1; e2222+=e2222Set[i];
                e1112Set[i]/=dataLicznik*0.1; e1112+=e1112Set[i];
                e1222Set[i]/=dataLicznik*0.1; e1222+=e1222Set[i];
            }
            e1111*=0.1; e1122*=0.1; e1212*=0.1; e2222*=0.1; e1112*=0.1; e1222*=0.1;
            //obliczenie bledow iloczynow elementow tensora odkształceń
            double dE1111=0, dE1122=0, dE1212=0, dE2222=0, dE1112=0, dE1222=0;
            for (int i=0;i<10;i++) {double epsilon=e1111-e1111Set[i]; dE1111+=epsilon*epsilon;} dE1111=getAvErrorFromSumEps(dE1111,90.0); //10*9 (n(n-1))
            for (int i=0;i<10;i++) {double epsilon=e1122-e1122Set[i]; dE1122+=epsilon*epsilon;} dE1122=getAvErrorFromSumEps(dE1122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1212-e1212Set[i]; dE1212+=epsilon*epsilon;} dE1212=getAvErrorFromSumEps(dE1212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e2222-e2222Set[i]; dE2222+=epsilon*epsilon;} dE2222=getAvErrorFromSumEps(dE2222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1112-e1112Set[i]; dE1112+=epsilon*epsilon;} dE1112=getAvErrorFromSumEps(dE1112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1222-e1222Set[i]; dE1222+=epsilon*epsilon;} dE1222=getAvErrorFromSumEps(dE1222,90.0);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... "); fflush(stdout);
            double s1111=e1111*avVolume, dS1111=fabs(e1111*dAvVolume)+fabs(dE1111*avVolume),
                   s1122=e1122*avVolume, dS1122=fabs(e1122*dAvVolume)+fabs(dE1122*avVolume),
                   s1212=e1212*avVolume, dS1212=fabs(e1212*dAvVolume)+fabs(dE1212*avVolume),
                   s2222=e2222*avVolume, dS2222=fabs(e2222*dAvVolume)+fabs(dE2222*avVolume),
                   s1112=e1112*avVolume, dS1112=fabs(e1112*dAvVolume)+fabs(dE1112*avVolume),
                   s1222=e1222*avVolume, dS1222=fabs(e1222*dAvVolume)+fabs(dE1222*avVolume),

                   nu2211_1111=-s1122/s1111, dNu2211_1111=fabs(dS1122/s1111)+fabs(dS1111*s1122/s1111/s1111),
                   nu1122_2222=-s1122/s2222, dNu1122_2222=fabs(dS1122/s2222)+fabs(dS2222*s1122/s2222/s2222),
                   avNu=(nu2211_1111+nu1122_2222)/2.0, dAvNu=(dNu2211_1111+dNu1122_2222)/2.0,

                   //\lambdaReduced=\lambda*\sigma^2/(kT)
                   l11=multimerS*multimerS/(8.0*(s1111+s1122)), dL11=multimerS*multimerS*(fabs(dS1111)+fabs(dS1122))/(8.0*fabs(s1111+s1122)*fabs(s1111+s1122)),
                   l12=multimerS*multimerS/(8.0*(s2222+s1122)), dL12=multimerS*multimerS*(fabs(dS2222)+fabs(dS1122))/(8.0*fabs(s2222+s1122)*fabs(s2222+s1122)),

                   B1=4.0*l11, dB1=4.0*dL11,
                   B2=4.0*l12, dB2=4.0*dL12,
                   avB=(B1+B2)/2.0, dAvB=(dB1+dB2)/2.0,

                   my1=(B1-B1*avNu)/(1.0+avNu), dMy1=fabs(dB1*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-B1/(1.0+avNu)-(B1-B1*avNu)/(1.0+avNu)/(1.0+avNu))),
                   my2=(B2-B2*avNu)/(1.0+avNu), dMy2=fabs(dB2*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-B2/(1.0+avNu)-(B2-B2*avNu)/(1.0+avNu)/(1.0+avNu))),
                   avMy=(my1+my2)/2.0, dAvMy=(dMy1+dMy2)/2.0,

                   avE=4.0*avB*avMy/(avB+avMy), dAvE=4.0*(fabs(avB*avB*dAvMy)+fabs(dAvB*avMy*avMy))/fabs(avB+avMy)/fabs(avB+avMy);
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla 1 czastki
            printf("Creation of 1-particle orientation file for Origin... "); fflush(stdout);
            double componentCounter=0, averageCos6PhiOne=0, ODFMaxOne=0;
            fileOrientations=fopen(bufferOrientations,"rt");
            double ODF_1P[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_1P[i]=0; //Orientational Distribution Function (1 Particle)
            int licznik=1,character;
            while ((character=fgetc(fileOrientations))!=EOF) {
                if (character==',') licznik++;
                if (licznik==3) {
                    licznik=0;
                    char data[50]=""; int actIndex=0;
                    while (true) {
                        character=fgetc(fileOrientations);
                        if (character!='}' && actIndex<50) data[actIndex++]=character;
                        else {
                            double angle = normalizeAngle(strtod(data,NULL)+C);
                            averageCos6PhiOne+=cos(6.0*angle); componentCounter++;
                            int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                            ODF_1P[index]++;
                            break;
                        }
                    }
                }
            }
            fclose(fileOrientations);
            double dPhi=2.0*C/((double)(ODFLength-1.0));
            double suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_1P[i]; for (int i=0;i<ODFLength;i++) ODF_1P[i]/=suma*dPhi;
            averageCos6PhiOne/=componentCounter; for (int i=0;i<ODFLength;i++) if (ODFMaxOne<ODF_1P[i]) ODFMaxOne=ODF_1P[i];
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla wszystkich czastek
            printf("Creation of ALL-particle orientation file for Origin... "); fflush(stdout);
            componentCounter=0; double averageCos6PhiAll=0, ODFMaxAll=0, averagePhiAll=0;
            fileAllOrientations = fopen(allOrientationsFileName,"rt");
            double ODF_AllP[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]=0; //Orientational Distribution Function (All Particles)
            char data[50]=""; int actIndex=0;
            while ((character=fgetc(fileAllOrientations))!=EOF) {
                if (character!=',' && character!='\n') data[actIndex++]=character;
                else {
                    data[actIndex++]=' '; actIndex=0;
                    double angle = normalizeAngle(strtod(data,NULL)+C);
                    averagePhiAll+=angle;
                    averageCos6PhiAll+=cos(6.0*angle); componentCounter++;
                    int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                    ODF_AllP[index]++;
                }
            }
            fclose(fileAllOrientations);
            suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_AllP[i]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]/=suma*dPhi;
            averagePhiAll/=componentCounter; averageCos6PhiAll/=componentCounter;
            int maxODFAllIndex;
            for (int i=0;i<ODFLength;i++) if (ODFMaxAll<ODF_AllP[i]) {
                ODFMaxAll=ODF_AllP[i];
                maxODFAllIndex=i;
            }
            printf("done\n");

            //analiza konfiguracji przejsciowych (dPhi z konfiguracji na konfiguracje)
            double avAbsDPhi=0;
            if (saveConfigurations) {
                printf("Transient configurations analysis... "); fflush(stdout);
                fileSavedConfigurations = fopen(bufferSavedConfigurations,"rt");
                double prevCfg[activeN][3]; bool undefinedPrevCfg=true;
                int CfgQuantity=0,CfgFileQuantity=0;
                while ((character=fgetc(fileSavedConfigurations))!=EOF) {
                    if (character=='n') {//text: newCfgFile (after merging) [it's NOT \n -> it's n, first char of 'new']
                        undefinedPrevCfg=true;
                        while (fgetc(fileSavedConfigurations)!='\n');
                        continue;
                    } else CfgQuantity++;
                    while (fgetc(fileSavedConfigurations)!='[');
                    for (int i=0;i<activeN;i++) {
                        for (int j=0;j<3;j++) {
                            strcpy(data,""); actIndex=0;
                            while ((character=fgetc(fileSavedConfigurations))!=',' && character!=']') data[actIndex++]=character;
                            data[actIndex++]=' ';
                            if (j==2) {
                                for (int k=0;k<3;k++) fgetc(fileSavedConfigurations); //skip ",m["
                                if (!undefinedPrevCfg) {
                                    avAbsDPhi+=fabs(strtod(data,NULL)-prevCfg[i][j]);
                                }
                            }
                            prevCfg[i][j]=strtod(data,NULL);
                        }
                    }
                    if (undefinedPrevCfg) {undefinedPrevCfg=false; CfgFileQuantity++;}
                }
                fclose(fileSavedConfigurations);
                avAbsDPhi/=(double)activeN*(CfgQuantity-CfgFileQuantity);
                printf("done\n");
            }

            long timeEndMath=time(0);




/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... "); fflush(stdout);
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) {
                fileConfigurations = fopen(bufferConfigurations,"w");
                fileConfigurationsList = fopen(configurationsListFileName,"a");
            }
            fileOrientationsResults = fopen(bufferOrientationsResults,"w");
            fileAllOrientationsResults = fopen(allOrientationsResultsFileName,"w");

            if (saveConfigurations) {
                fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%ld\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avVolume,avBoxMatrix[0],avBoxMatrix[1],avBoxMatrix[2],avRho,avPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,s1112,dS1112,s1222,dS1222,avNu,dAvNu,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,savedConfigurationsInt,avAbsDPhi);
                fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%ld\t%.17E\n",pressureReduced,avPacFrac,avNu,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,avB,avMy,avE,savedConfigurationsInt,avAbsDPhi);
            } else {
                fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avVolume,avBoxMatrix[0],avBoxMatrix[1],avBoxMatrix[2],avRho,avPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,s1112,dS1112,s1222,dS1222,avNu,dAvNu,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll);
                fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",pressureReduced,avPacFrac,avNu,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,avB,avMy,avE);
            }

            if (!onlyMath[0]) {
                rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                fprintf(fileConfigurationsList,"multimers[x_,y_,kI_]:={");
                for (int i=0;i<activeN;i++) {
                    fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                    fprintf(fileConfigurationsList,"m[%.12E+x,%.12E+y,%.12E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                }
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                fprintf(fileConfigurationsList,"{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[If[x==0 && y==0,0.4,0]],Green,Disk[{%.12E,%.12E},%.12E]}};\nconfigurationsList=Append[configurationsList,g[%.12E,multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                        boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,pacFrac,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
            }

            for (int i=0;i<ODFLength;i++) {
                fprintf(fileOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_1P[i]);
                fprintf(fileAllOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_AllP[i]);
            }

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) {
                fclose(fileConfigurations); fclose(fileConfigurationsList);
            }
            fclose(fileOrientationsResults); fclose(fileAllOrientationsResults);
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        getNextArgument(arg,true);
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        args[4]=0;
        loadedConfiguration=0;
        generatorStartPoint=0;
    }
    printf("\nTime for equilibrations: %ldsec ,  time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n");
}
