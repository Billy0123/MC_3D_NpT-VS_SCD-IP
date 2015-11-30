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
    onlyMath[2]={0,0};
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,intervalOrientations,savedConfigurationsInt;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg,deltaR,deltaPhi,deltaV=0.1, //*sigma(=1)
       multimerS,multimerD,randomStartStep[2],normalizingDenominator,
       neighRadius,neighRadius2,neighRadiusMod,neighSafeDistance,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][2],pi=3.1415926535898;
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
    int neighbours[50], neighCounter;
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

int initPositions (particle *particles, double boxMatrix[2][2], double matrixOfParticlesSize[2], double rowShift, double initPeriodicImageMinDistance, double startAngleInRows[2]) {
    double mod=sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]), interval[2], actualPosition[2]={0,0};
    for (int i=0;i<2;i++) interval[i]=boxMatrix[i][i]/matrixOfParticlesSize[i]/mod;
    int rowCounter=0;
    for (int i=0;i<N;i++) {
        for (int j=0;j<2;j++) particles[i].r[j]=actualPosition[j];
        particles[i].normR[0]=(-boxMatrix[1][1]*particles[i].r[0]+boxMatrix[1][0]*particles[i].r[1])/normalizingDenominator;
        particles[i].normR[1]=(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/normalizingDenominator;
        particles[i].phi=startAngleInRows[rowCounter%2];
        actualPosition[0]+=interval[0];
        if ((rowCounter%2)*rowShift-actualPosition[0]+boxMatrix[0][0]<initPeriodicImageMinDistance) {  //sprawdzenie czy kolejna wstawiona cząstka nie przekryje swojego obrazu periodycznego
            actualPosition[0]=(++rowCounter%2)*rowShift;
            actualPosition[1]+=interval[1];
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
                        particles[k].phi+=0.000001;
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
            if (dr<minDistance) energy=1; //analyticMethodForEvenHCM 1/6
            else {
                double gamma=atan(ry/rx),
                       aAngle=particles[index].phi-gamma,
                       bAngle=particles[particles[index].neighbours[i]].phi-gamma;
                if (multimerN%2!=0) {
                    if (rx>0) bAngle-=C;
                    else aAngle-=C;
                }
                aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle); //analyticMethodForEvenHCM 2/6
                energy=checkOverlaps(dr,aAngle,bAngle);
                //energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
            } //analyticMethodForEvenHCM 3/6
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
    particles[index].normR[0]=(-boxMatrix[1][1]*particles[index].r[0]+boxMatrix[1][0]*particles[index].r[1])/normalizingDenominator;
    particles[index].normR[1]=(boxMatrix[1][0]*particles[index].r[0]-boxMatrix[0][0]*particles[index].r[1])/normalizingDenominator;
    particles[index].phi+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaPhi;
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
    double lnNewBoxMatrix[2][2], newBoxMatrix[2][2];
    if (pressure>pressureRealOfNotFluid) {  //dozwolone zmiany ksztaltu pudla (faza stala)
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) lnNewBoxMatrix[i][j]=log(boxMatrix[i][j]);
        lnNewBoxMatrix[0][0]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        lnNewBoxMatrix[1][1]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        if (boxMatrix[1][0]==0) lnNewBoxMatrix[1][0]=-20.0+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV; //log(0) -> -\infty; E^~0 za duze jak na start
        else lnNewBoxMatrix[1][0]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        lnNewBoxMatrix[0][1]=lnNewBoxMatrix[1][0];
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) newBoxMatrix[i][j]=exp(lnNewBoxMatrix[i][j]);
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        lnNewBoxMatrix[0][0]=log(boxMatrix[0][0])+(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaV;
        newBoxMatrix[0][0]=exp(lnNewBoxMatrix[0][0]);
        newBoxMatrix[1][1]=boxMatrix[1][1]/boxMatrix[0][0]*newBoxMatrix[0][0];
        newBoxMatrix[0][1]=boxMatrix[0][1]; newBoxMatrix[1][0]=boxMatrix[1][0];
    }
    double newVolume=fabs(-pow(newBoxMatrix[1][0],2)+newBoxMatrix[0][0]*newBoxMatrix[1][1]);

    double bufferR[activeN][2];
    for (int i=0;i<activeN;i++) {
        bufferR[i][0]=newBoxMatrix[0][0]*particles[i].normR[0]+newBoxMatrix[1][0]*particles[i].normR[1];
        bufferR[i][1]=newBoxMatrix[0][1]*particles[i].normR[0]+newBoxMatrix[1][1]*particles[i].normR[1];
    }

    //matrix 11/16
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
                if (dr<minDistance) energy=1; //analyticMethodForEvenHCM 4/6
                else {
                    double gamma=atan(ry/rx),
                           aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                           bAngle=particles[i].phi-gamma;
                    if (multimerN%2!=0) {
                        if (rx>0) bAngle-=C;
                        else aAngle-=C;
                    }
                    aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle); //analyticMethodForEvenHCM 5/6
                    energy=checkOverlaps(dr,aAngle,bAngle);
                    //energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
                } //analyticMethodForEvenHCM 6/6
                if (energy==1) {
                    result=0;
                    i=activeN; break;
                }
            }
        }
    }

    if (result) {
        double arg=-(pressure*(newVolume-(*volume))-(((double)N+1)*log(newVolume/(*volume))+log((newBoxMatrix[0][0]+newBoxMatrix[1][1])/(boxMatrix[0][0]+boxMatrix[1][1]))));
        if (MTGenerate(randomStartStep)%1000000/1000000.0>exp(arg)) result=0;
        if (result) {
            *volume=newVolume;
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
            for (int i=0;i<activeN;i++) for (int j=0;j<2;j++) particles[i].r[j]=bufferR[i][j];
            normalizingDenominator=pow(boxMatrix[1][0],2)-boxMatrix[0][0]*boxMatrix[1][1];
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
        iterationTable[fileIterateIterationsNumber++][1]=strtod(pEnd,NULL);
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

double getNextArgument (double prevArg, bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            prevArg=growing?iterationTable[actIteration][1]:iterationTable[fileIterateIterationsNumber-1-actIteration][1];
            if (growing) startMinPacFrac=iterationTable[actIteration][0];
            else startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg*10000)>=round(intervalMin[i]*10000) && round(prevArg*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg*10000)+round(intervalDelta[i]*10000);
                prevArg=newArg/10000.0;
                break;
            }
        if (round(prevArg*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg*10000)>round(intervalMin[i]*10000) && round(prevArg*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg*10000)-round(intervalDelta[i]*10000);
                prevArg=newArg/10000.0;
                break;
            }
        if (round(prevArg*10000)<round(minArg*10000)) growing=-1;
    }
    return prevArg;
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
        while(fgets(config,300,fileConfig)!=NULL) {
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
                case 4:multimerS=strtod(data,NULL);break;
                case 5:multimerD=strtod(data,NULL);break;
                case 6:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 7:growing=strtol(data,NULL,10);break;
                case 8:loadedConfiguration=strtol(data,NULL,10);break;
                case 9:loadedArg=strtod(data,NULL);break;
                case 10:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 11:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 12:loadedSetGenerator=strtol(data,NULL,10);break;
                case 13:iterationsNumber=strtol(data,NULL,10);break;
                case 14:countCollidingPairs=strtol(data,NULL,10);break;
                case 15:intervalSampling=strtol(data,NULL,10);break;
                case 16:intervalOutput=strtol(data,NULL,10);break;
                case 17:saveConfigurations=strtol(data,NULL,10);break;
                case 18:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 19:ODFLength=strtol(data,NULL,10);break;
                case 20:OCFMode=strtol(data,NULL,10);break;
                case 21:neighRadiusMod=strtod(data,NULL);break;
                case 22:intervalOrientations=strtol(data,NULL,10);break;
                case 23:skipFirstIteration=strtol(data,NULL,10);break;
                case 24:useSpecificDirectory=strtol(data,NULL,10);break;
                case 25:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 26:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 27:intervalResults=strtol(data,NULL,10);break;
                case 28:maxDeltaR=strtod(data,NULL);break;
                case 29:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 30:desiredAcceptanceRatioV=strtod(data,NULL);break;
                case 31:useFileToIterate=strtol(data,NULL,10);break;
                case 32:startMinPacFrac=strtod(data,NULL);break;
                case 33:startMaxPacFrac=strtod(data,NULL);break;
                case 34:minArg=strtod(data,NULL);break;
                case 35:maxArg=strtod(data,NULL);break;
                case 36:multiplyArgument=strtol(data,NULL,10);break;
                case 37:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-38)%3) {
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
                if (argumentsNumber==12) {
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
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==12) {
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
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==13) {
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
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==13) {
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
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==12) {
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
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        if (growing) {
            startMinPacFrac=iterationTable[0][0];
            startArg=iterationTable[0][1];
        } else {
            startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1][0];
            startArg=iterationTable[fileIterateIterationsNumber-1][1];
        }
    } else startArg=growing?minArg:maxArg;
    pressureRealOfNotFluid/=(multimerS*multimerS);
    deltaR=maxDeltaR*multimerS; deltaV*=multimerS;
    for (int i=0;i<pointNumber;i++) startArg=getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg;
    activeN=N-gaps;

    //sprawdzenie czy uzyskana konfiguracja jest obsługiwana
    if (multimerN!=6 && multimerN!=5) {
        printf("ERROR: Not supported multimerN: %d.\n",multimerN);
        return 0;
    }
    if (N%56!=0 && N%780!=0) {
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
    if (multimerN==6) VcpPerParticle=minDistance*minDistance*sqrt(3)/2.0;  //dla heksamerow o dowolnym d/\sigma
    else if (multimerN==5) VcpPerParticle=5.0936*multimerD*multimerD;  //dla pentamerów o d/\sigma=1
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
    long arg5=0; double arg1, arg2, arg3, arg4, arg6, arg7, arg8, arg9, arg10;

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileOrientations, *fileOrientatCorrelFun, *fileConfigurationsList, *fileAllResults, *fileAllOrientations, *fileOrientationsResults, *fileAllOrientationsResults;
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        if (saveConfigurations) fprintf(fileResults,"Cycles\tPressureReduced\tVolume\tBoxMatrix[0][0]\tBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tRho\tV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tavNu\tdAvNu\tavB\tdAvB\tavMy\tdAvMy\tavE\tdAvE\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileResults,"Cycles\tPressureReduced\tVolume\tBoxMatrix[0][0]\tBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tRho\tV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tavNu\tdAvNu\tavB\tdAvB\tavMy\tdAvMy\tavE\tdAvE\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        if (saveConfigurations) fprintf(fileExcelResults,"PressureReduced\tV/V_cp\tavNu\tODFMax_All\t<cos(6Phi)>_All\tavB\tavMy\tavE\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileExcelResults,"PressureReduced\tV/V_cp\tavNu\tODFMax_All\t<cos(6Phi)>_All\tavB\tavMy\tavE\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg=startArg, oldVolume, oldBoxMatrix[2][2];
    while (growing>=0) {
        double pressureReduced=arg, pressureReal=pressureReduced/multimerS/multimerS, //pressureReduced=\tau*\sigma^2/(kT), kT=1
               volume=(startArg==arg)?(growing?N/(1.0/VcpPerParticle/startMinPacFrac):N/(1.0/VcpPerParticle/startMaxPacFrac)):oldVolume, rho=N/volume, pacFrac=1.0/VcpPerParticle/rho,
               boxMatrix[2][2],unitCellAtCP[2],initRowShift,initPeriodicImageMinDistance,startAngleInRows[2],matrixOfParticlesSize[2];
        if (multimerN==6) {
            unitCellAtCP[0]=minDistance; unitCellAtCP[1]=sqrt(3)/2.0*minDistance;
            initRowShift=unitCellAtCP[0]*sqrt(pacFrac)/2.0; initPeriodicImageMinDistance=minDistance;
            startAngleInRows[0]=absoluteMinimum2; startAngleInRows[1]=absoluteMinimum2;
        } else if (multimerN==5) {
            unitCellAtCP[0]=2.40487*multimerD; unitCellAtCP[1]=2.11804*multimerD;
            initRowShift=1.39176*multimerD; initPeriodicImageMinDistance=getMinimalDistanceAnalyticalMethodForOddHCM(0,0);
            startAngleInRows[0]=0; startAngleInRows[1]=C;
        }
        if (N%56==0) {matrixOfParticlesSize[0]=7; matrixOfParticlesSize[1]=8;}
        else if (N%780==0) {matrixOfParticlesSize[0]=26; matrixOfParticlesSize[1]=30;}
        double NLinearMod = sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]);
        if (startArg==arg) {
            for (int i=0;i<2;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]*sqrt(pacFrac)*NLinearMod;
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0;
        } else for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
        normalizingDenominator=pow(boxMatrix[1][0],2)-boxMatrix[0][0]*boxMatrix[1][1];

        if (!onlyMath[0]) {
            if (arg==startArg && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.4f (StartDen: %.4f, startPacFrac: %.4f), mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,startArg,rho,pacFrac,multimerN,multimerS,multimerD);
                if (!initPositions(particles,boxMatrix,matrixOfParticlesSize,initRowShift,initPeriodicImageMinDistance,startAngleInRows)) return 0;
                adjustAngles(particles,boxMatrix);
                updateNeighbourList(particles,boxMatrix); //matrix(comment) 12/16
            } else if (loadedConfiguration) {
                char configurations[400+110*activeN];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                int numerLinii=1;
                while(fgets(configurations,400+110*activeN,fileCTL)!=NULL && numerLinii++<=3) sscanf(configurations,"%c",configurations);
                fclose(fileCTL);

                char *pEnd;
                arg1=strtod(configurations,&pEnd);  //RandStart
                arg2=strtod(pEnd,&pEnd);            //RandStep
                arg3=strtod(pEnd,&pEnd);            //Density
                arg4=strtod(pEnd,&pEnd);            //PressureReduced
                arg5=strtol(pEnd,&pEnd,10);         //Cycles
                arg6=strtod(pEnd,&pEnd);            //boxMatrix[0][0]
                arg7=strtod(pEnd,&pEnd);            //boxMatrix[1][1]
                arg8=strtod(pEnd,&pEnd);            //boxMatrix[0][1]=boxMatrix[1][0] (symetryczna macierz)
                arg9=strtod(pEnd,&pEnd);            //deltaR
                arg10=strtod(pEnd,NULL);            //deltaV

                boxMatrix[0][0]=arg6; boxMatrix[1][1]=arg7; boxMatrix[1][0]=arg8; boxMatrix[0][1]=arg8;
                deltaR=arg9; deltaPhi=deltaR*2.0*sin(C); deltaV=arg10;
                arg=arg4; pressureReduced=arg; pressureReal=pressureReduced/multimerS/multimerS;
                rho=arg3; pacFrac=1.0/VcpPerParticle/rho; volume=N/rho;

                normalizingDenominator=pow(boxMatrix[1][0],2)-boxMatrix[0][0]*boxMatrix[1][1];
                int actIndex=0;
                while (configurations[actIndex]!='{') actIndex++;
                actIndex+=3;
                for (int i=0;i<activeN;i++) {
                    for (int j=0;j<3;j++) {
                        char coordinate[50];
                        int licznik=0;
                        while (configurations[actIndex]!=',') coordinate[licznik++]=configurations[actIndex++];
                        if (j<2) {
                            actIndex++;
                            particles[i].r[j]=strtod(coordinate,NULL);
                        } else {
                            //actIndex+=22;         //dla mniejszej dokladnosci (%.5f) - dla kompatybilnosci ze starymi plikami
                            actIndex+=36;
                            particles[i].phi=strtod(coordinate,NULL);
                        }
                    }
                    particles[i].normR[0]=(-boxMatrix[1][1]*particles[i].r[0]+boxMatrix[1][0]*particles[i].r[1])/normalizingDenominator;
                    particles[i].normR[1]=(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/normalizingDenominator;
                }
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.4f (startDen: %.4f, startPacFrac: %.4f), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4f, DeltaV: %.4f\n",N,gaps,growing,arg4,arg3,pacFrac,arg1,arg2,arg5,arg9,arg10);
                //for (int i=0;i<2;i++) for (int j=0;j<2;j++) printf("boxMatrix[%d][%d]=%.12f\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<N;i++) printf("%.12f,  %.12f,  %.12f\n",particles[i].r[0],particles[i].r[1],particles[i].phi);return 0;
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
                        InitMT((unsigned int)arg1);
                        randomStartStep[0]=arg1;
                    } else {
                        randomStartStep[0]=InitRandomMT();
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                    }
                    randomStartStep[1]=0;
                    if (loadedSetGenerator) {
                        printf("Setting p-random number generator to last position from file...\n");
                        for (double i=0;i<arg2;i++) MTGenerate(randomStartStep);
                    } else printf("Setting p-random number generator to last position from file - DISABLED\n");
                } else {
                    printf("Setting start position of p-random number generator to actual CPU time...\n");
                    randomStartStep[0]=InitRandomMT();
                    randomStartStep[1]=0;
                }
                printf("Start of equilibration at reduced pressure: %.4f (startDen: %.4f, startPacFrac: %.4f)... (%ld cycles)\n",pressureReduced,rho,pacFrac,cyclesOfEquilibration);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.4f, mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,pressureReduced,multimerN,multimerS,multimerD);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                fullCycle=activeN+volumeMoveChance,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberV=0, displacedNumberV=0;
            double nStep=fullCycle*((double)cyclesOfEquilibration+(double)cyclesOfMeasurement+10.0),
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
                        possibleDistance+=sqrt(0.5*deltaR*deltaR)*((double)intervalSampling)*acceptanceRatioR*5.0;  //ostatnie *5.0 - dla bezpieczenstwa; sqrt(0.5*deltaR*0.5*deltaR+0.5*deltaR*0.5*deltaR)=sqrt(0.5*deltaR*deltaR)
                        if (possibleDistance>=neighSafeDistance) {
                            updateNeighbourList(particles,boxMatrix);
                            possibleDistance=0;
                        }  //matrix(comment) 16/16

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
                                            if (energy==1) printf("colliding: distance- %.12f, alpha: %.12f, beta: %.12f, analytic: %.12f, i: %d, j: %d\n",dr,aAngle,bAngle,getMinimalDistanceAnalyticalMethodForEvenHCM(normalizeAngle(aAngle+C),normalizeAngle(bAngle+C)),i,j);
                                        }
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+arg5),collidingPairs);
                            else printf("Cycle: %ld\n",(cycle+arg5));
                            printf("   AccRatR: %.4f, dR: %.4f, AccRatV: %.4f, dV: %.4f\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                            printf("   Dens: %.4f, V/V_cp: %.4f, PressRed: %.4f\n",rho,pacFrac,pressureReduced);
                            printf("   box00: %.8f, box11: %.8f, box10(01): %.8f\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
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
                                        else {
                                            double gamma=atan(ry/rx),
                                                    aAngle=particles[j].phi-gamma,
                                                    bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (rx>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);  //analyticCheck 2/3
                                            energy=checkOverlaps(dr,aAngle,bAngle);
                                            //energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
                                        }  //analyticCheck 3/3
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0)
                                fprintf(fileAllResults,"%ld\t%.17f\t%.17f\t%.17f\t%.17f\t%.17f\t%.17f\t\n",(cycle+arg5),volume,boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0],rho,pacFrac);

                            if (cycle%intervalOrientations==0) {
                                /////skan po konkretnej cząstce (np. dla 224: 14*(16/2)-(14/2)=105, etc.) - w srodku by nie skakala na granicy pudla periodycznego; bezposrednie uzycie w Mathematice (format tablicy)
                                if (cycle-cyclesOfEquilibration>=cyclesOfMeasurement)
                                    fprintf(fileOrientations,"{%.12f,%.12f,%.12f}}",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                else fprintf(fileOrientations,"{%.12f,%.12f,%.12f},",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                /////skan po wszystkich cząstkach
                                for (int i=0;i<activeN-1;i++) fprintf(fileAllOrientations,"%.12f,",particles[i].phi);
                                fprintf(fileAllOrientations,"%.12f\n",particles[activeN-1].phi);
                                /////OCF
                                if (OCFMode) {
                                    char bufferText[1000000]="", bufferAngle[20];
                                    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                        if (i<particles[i].neighbours[j]) {
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

                                            strncat(bufferText,"{",2); sprintf(bufferAngle,"%.12f",aAngle); strncat(bufferText,bufferAngle,20);
                                            strncat(bufferText,",",2); sprintf(bufferAngle,"%.12f",bAngle); strncat(bufferText,bufferAngle,20);
                                            strncat(bufferText,"},",3);
                                        }
                                    }
                                    if (cycle-cyclesOfEquilibration>=cyclesOfMeasurement) bufferText[strlen(bufferText)-1]='}';
                                    fprintf(fileOrientatCorrelFun,"%s",bufferText);
                                }
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.12f\t%.12f\t%.12f\t{",(cycle+arg5),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++)
                                    fprintf(fileSavedConfigurations,"m[%.17f,%.17f,%.17f],",particles[i].r[0],particles[i].r[1],particles[i].phi);
                                fprintf(fileSavedConfigurations,"m[%.17f,%.17f,%.17f]}",particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].phi);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+arg5),collidingPairs);
                                else printf("Cycle: %ld\n",(cycle+arg5));
                                printf("   AccRatR: %.4f, dR: %.4f, AccRatV: %.4f, dV: %.4f\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                                printf("   Dens: %.4f, V/V_cp: %.4f, PressRed: %.4f\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8f, box11: %.8f, box10(01): %.8f\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12f\tV/V_cp: %.12f\tPressureRed: %.12f\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+arg5),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17f %.17f %ld %.17f %.17f %.17f %.17f %.17f {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+arg5),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                                for (int i=0;i<activeN;i++)
                                    fprintf(fileConfigurations,"m[%.17f,%.17f,%.17f,%.12f,%.12f,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12f,%.12f},{%.12f,%.12f},{%.12f,%.12f}}]},{Opacity[0.2],Green,Disk[{%.12f,%.12f},%.12f]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12f, boxMatrix[1][1]=%.12f, boxMatrix[1][0]=boxMatrix[0][1]=%.12f",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                                fclose(fileConfigurations);
                            }
                        } else {
                            if (acceptanceRatioR>desiredAcceptanceRatioR) deltaR*=1.05; else deltaR*=0.95;
                            if (deltaR>maxDeltaR) deltaR=maxDeltaR;
                            deltaPhi=deltaR*2.0*sin(C);
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; else deltaV*=0.95;
                        }
                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberV=0; displacedNumberV=0;
                    }
                    if (cycle-cyclesOfEquilibration>=cyclesOfMeasurement) iStep=nStep;
                }
            }
            fclose(fileAllResults);
            fclose(fileOrientations); fclose(fileAllOrientations);
            if (saveConfigurations) fclose(fileSavedConfigurations); if (OCFMode) fclose(fileOrientatCorrelFun);
            if (timeEquilibration==0) timeEquilibration=time(0);
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            printf("Start of calculation of results...\n");

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... ");
            double avVolume=0, avBoxMatrix[3]={0,0,0}, avRho=0, avPacFrac=0;
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[activeN*17];
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);

            long dataLicznik=0;
            while(fgets(linia,300,fileAllResults)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++; if (actIndex>=300) continue;
                int dataIndex=0; double dataD[6]; while (dataIndex<6) {
                    char data[50];
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    dataD[dataIndex++]=strtod(data,NULL);
                }
                if (dataIndex<10) {
                    avVolume+=dataD[0]; avBoxMatrix[0]+=dataD[1];
                    avBoxMatrix[1]+=dataD[2]; avBoxMatrix[2]+=dataD[3];
                    avRho+=dataD[4]; avPacFrac+=dataD[5];
                    dataLicznik++;
                }
            }
            fclose(fileAllResults);
            avVolume/=(double)dataLicznik; avRho/=(double)dataLicznik; avPacFrac/=(double)dataLicznik;
            for (int i=0;i<3;i++) avBoxMatrix[i]/=(double)dataLicznik;
            printf("done\n");

            //obliczenie bledu objetosci (dV)
            printf("Calculation of dV... ");
            double dAvVolume = 0;
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            while(fgets(linia,300,fileAllResults)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++;
                if (actIndex>=300) continue;
                char data[50];
                int licznik=0;
                while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) continue;
                double epsilon=avVolume-strtod(data,NULL);
                dAvVolume+=epsilon*epsilon;
            }
            fclose(fileAllResults);
            double denominator = (double)dataLicznik*((double)dataLicznik-1.0);
            dAvVolume=getAvErrorFromSumEps(dAvVolume,denominator);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... ");
            double e1111=0, e1122=0, e1212=0, e2222=0,
                   HxyHyx=avBoxMatrix[2]*avBoxMatrix[2],
                   HxxHyy=avBoxMatrix[0]*avBoxMatrix[1],
                   HxxHxy=avBoxMatrix[0]*avBoxMatrix[2],
                   Hxx2=avBoxMatrix[0]*avBoxMatrix[0],
                   Hyy2=avBoxMatrix[1]*avBoxMatrix[1],
                   mod0=HxyHyx-HxxHyy,
                   mod1=1.0/(2.0*mod0*mod0);
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            while(fgets(linia,300,fileAllResults)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++;
                if (actIndex>=300) continue;
                double h11,h22,h12;
                int dataIndex=0; while (dataIndex<4) {
                    char data[50];
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    switch (dataIndex++) {
                        case 1: h11=strtod(data,NULL);break;
                        case 2: h22=strtod(data,NULL);break;
                        case 3: h12=strtod(data,NULL);break;
                    }
                }
                if (dataIndex<10) {
                    double hxyhyx=h12*h12,
                           hxxPhyy=h11+h22,
                           hxx2=h11*h11,
                           hyy2=h22*h22,

                           e11=mod1*(HxyHyx*(hxyhyx-HxyHyx+hyy2)-2.0*(h11*avBoxMatrix[2]*h12-avBoxMatrix[0]*HxyHyx+avBoxMatrix[2]*h12*h22)*avBoxMatrix[1]+(hxx2-Hxx2+hxyhyx)*Hyy2),
                           e22=mod1*(HxyHyx*(hxx2+hxyhyx-HxyHyx)-2.0*(HxxHxy*h12*hxxPhyy-HxxHyy*HxyHyx)+Hxx2*(hxyhyx+hyy2-Hyy2)),
                           e12=mod1*(-HxxHxy*(hxyhyx+hyy2)+HxxHyy*h12*hxxPhyy+avBoxMatrix[2]*(h12*avBoxMatrix[2]*hxxPhyy-(hxx2+hxyhyx)*avBoxMatrix[1]));

                    e1111+=e11*e11;
                    e1122+=e11*e22;
                    e1212+=e12*e12;
                    e2222+=e22*e22;
                }
            }
            fclose(fileAllResults);
            e1111/=(double)dataLicznik; e1122/=(double)dataLicznik;
            e1212/=(double)dataLicznik; e2222/=(double)dataLicznik;
            printf("done\n");

            //obliczenie bledow iloczynow elementow tensora odksztalcen
            printf("Calculation of errors of average values of products of strain tensor's elements... ");
            double dE1111=0, dE1122=0, dE1212=0, dE2222=0;
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            while(fgets(linia,300,fileAllResults)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++;
                if (actIndex>=300) continue;
                double h11,h22,h12;
                int dataIndex=0; while (dataIndex<4) {
                    char data[50];
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    switch (dataIndex++) {
                        case 1: h11=strtod(data,NULL);break;
                        case 2: h22=strtod(data,NULL);break;
                        case 3: h12=strtod(data,NULL);break;
                    }
                }
                if (dataIndex<10) {
                    double hxyhyx=h12*h12,
                           hxxPhyy=h11+h22,
                           hxx2=h11*h11,
                           hyy2=h22*h22,

                           e11=mod1*(HxyHyx*(hxyhyx-HxyHyx+hyy2)-2.0*(h11*avBoxMatrix[2]*h12-avBoxMatrix[0]*HxyHyx+avBoxMatrix[2]*h12*h22)*avBoxMatrix[1]+(hxx2-Hxx2+hxyhyx)*Hyy2),
                           e22=mod1*(HxyHyx*(hxx2+hxyhyx-HxyHyx)-2.0*(HxxHxy*h12*hxxPhyy-HxxHyy*HxyHyx)+Hxx2*(hxyhyx+hyy2-Hyy2)),
                           e12=mod1*(-HxxHxy*(hxyhyx+hyy2)+HxxHyy*h12*hxxPhyy+avBoxMatrix[2]*(h12*avBoxMatrix[2]*hxxPhyy-(hxx2+hxyhyx)*avBoxMatrix[1]));

                    double epsilon=e1111-e11*e11; dE1111+=epsilon*epsilon;
                    epsilon=e1122-e11*e22; dE1122+=epsilon*epsilon;
                    epsilon=e1212-e12*e12; dE1212+=epsilon*epsilon;
                    epsilon=e2222-e22*e22; dE2222+=epsilon*epsilon;
                }
            }
            fclose(fileAllResults);
            dE1111=getAvErrorFromSumEps(dE1111,denominator); dE1122=getAvErrorFromSumEps(dE1122,denominator);
            dE1212=getAvErrorFromSumEps(dE1212,denominator); dE2222=getAvErrorFromSumEps(dE2222,denominator);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... ");
            double s1111=e1111*avVolume, dS1111=fabs(e1111*dAvVolume)+fabs(dE1111*avVolume),
                   s1122=e1122*avVolume, dS1122=fabs(e1122*dAvVolume)+fabs(dE1122*avVolume),
                   s1212=e1212*avVolume, dS1212=fabs(e1212*dAvVolume)+fabs(dE1212*avVolume),
                   s2222=e2222*avVolume, dS2222=fabs(e2222*dAvVolume)+fabs(dE2222*avVolume),

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
            printf("Creation of 1-particle orientation file for Origin... ");
            double componentCounter=0, averageCos6PhiOne=0, ODFMaxOne=0;
            fileOrientations=fopen(bufferOrientations,"rt");
            double ODF_1P[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_1P[i]=0; //Orientational Distribution Function (1 Particle)
            int licznik=1;
            while (fgets(linia,2,fileOrientations)!=NULL) {
                sscanf(linia,"%c",linia);
                if (linia[0]==',') licznik++;
                if (licznik==3) {
                    licznik=0;
                    char data[50]; int actIndex=0;
                    while (true) {
                        fgets(linia,2,fileOrientations); sscanf(linia,"%c",linia);
                        if (linia[0]!='}') data[actIndex++]=linia[0];
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
            printf("Creation of ALL-particle orientation file for Origin... ");
            componentCounter=0; double averageCos6PhiAll=0, ODFMaxAll=0;
            fileAllOrientations = fopen(allOrientationsFileName,"rt");
            double ODF_AllP[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]=0; //Orientational Distribution Function (All Particles)
            while(fgets(linia,activeN*50,fileAllOrientations)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0, licznikN=0;
                while (actIndex<activeN*30+10 && licznikN++<activeN) {
                    char data[50]; int licznik=0;
                    while (linia[actIndex]!=',' && licznik<30) data[licznik++]=linia[actIndex++]; actIndex++;
                    double angle = normalizeAngle(strtod(data,NULL)+C);
                    averageCos6PhiAll+=cos(6.0*angle); componentCounter++;
                    int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                    ODF_AllP[index]++;
                }
            }
            fclose(fileAllOrientations);
            suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_AllP[i]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]/=suma*dPhi;
            averageCos6PhiAll/=componentCounter; for (int i=0;i<ODFLength;i++) if (ODFMaxAll<ODF_AllP[i]) ODFMaxAll=ODF_AllP[i];
            printf("done\n");

            //analiza konfiguracji przejsciowych
            double avAbsDPhi=0;
            if (saveConfigurations) {
                printf("Transient configurations analysis... ");
                fileSavedConfigurations = fopen(bufferSavedConfigurations,"rt");
                char configurations[activeN*70+100];
                double prevCfg[activeN][3]; bool undefinedPrevCfg=true;
                int CfgQuantity=0,CfgFileQuantity=0;
                while(fgets(configurations,activeN*70+100,fileSavedConfigurations)!=NULL) {//newCfgFile
                    int actIndex=0;
                    if (configurations[actIndex]=='n') {    //newCfgFile (after merging)
                        undefinedPrevCfg=true;
                        continue;
                    } else CfgQuantity++;
                    while (configurations[actIndex]!='{') actIndex++;
                    actIndex+=3;
                    for (int i=0;i<activeN;i++) {
                        for (int j=0;j<3;j++) {
                            char coordinate[50];
                            licznik=0;
                            while (configurations[actIndex]!=',' && configurations[actIndex]!=']') coordinate[licznik++]=configurations[actIndex++];
                            if (j<2) {  //polozenie
                                actIndex++;
                            } else {    //orientacja
                                actIndex+=4;
                                if (!undefinedPrevCfg) {
                                    avAbsDPhi+=fabs(strtod(coordinate,NULL)-prevCfg[i][j]);
                                }
                            }
                            prevCfg[i][j]=strtod(coordinate,NULL);
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

            printf("Saving data to files... ");
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
                fprintf(fileResults,"%ld\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%ld\t%.12f\n",(cycle+arg5),pressureReduced,avVolume,avBoxMatrix[0],avBoxMatrix[1],avBoxMatrix[2],avRho,avPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,avNu,dAvNu,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,savedConfigurationsInt,avAbsDPhi);
                fprintf(fileExcelResults,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%ld\t%.12f\n",pressureReduced,avPacFrac,avNu,ODFMaxAll,averageCos6PhiAll,avB,avMy,avE,savedConfigurationsInt,avAbsDPhi);
            } else {
                fprintf(fileResults,"%ld\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",(cycle+arg5),pressureReduced,avVolume,avBoxMatrix[0],avBoxMatrix[1],avBoxMatrix[2],avRho,avPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,avNu,dAvNu,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll);
                fprintf(fileExcelResults,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",pressureReduced,avPacFrac,avNu,ODFMaxAll,averageCos6PhiAll,avB,avMy,avE);
            }

            if (!onlyMath[0]) {
                rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12f\tV/V_cp: %.12f\tPressureRed: %.12f\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+arg5),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17f %.17f %ld %.17f %.17f %.17f %.17f %.17f {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+arg5),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                fprintf(fileConfigurationsList,"multimers[x_,y_,kI_]:={");
                for (int i=0;i<activeN;i++) {
                    fprintf(fileConfigurations,"m[%.17f,%.17f,%.17f,%.12f,%.12f,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                    fprintf(fileConfigurationsList,"m[%.12f+x,%.12f+y,%.12f,%.12f,%.12f,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                }
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12f,%.12f},{%.12f,%.12f},{%.12f,%.12f}}]},{Opacity[0.2],Green,Disk[{%.12f,%.12f},%.12f]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12f, boxMatrix[1][1]=%.12f, boxMatrix[1][0]=boxMatrix[0][1]=%.12f",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                fprintf(fileConfigurationsList,"{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12f,%.12f},{%.12f,%.12f},{%.12f,%.12f}}]},{Opacity[If[x==0 && y==0,0.4,0]],Green,Disk[{%.12f,%.12f},%.12f]}};\nconfigurationsList=Append[configurationsList,g[%.12f,multimers[%.12f,%.12f,kolorIndex=1],multimers[%.12f,%.12f,kolorIndex=1],multimers[%.12f,%.12f,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                        boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,pacFrac,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
            }

            for (int i=0;i<ODFLength;i++) {
                fprintf(fileOrientationsResults,"%.12f\t%.12f\n",-C+i*dPhi,ODF_1P[i]);
                fprintf(fileAllOrientationsResults,"%.12f\t%.12f\n",-C+i*dPhi,ODF_AllP[i]);
            }

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) {
                fclose(fileConfigurations); fclose(fileConfigurationsList);
            }
            fclose(fileOrientationsResults); fclose(fileAllOrientationsResults);
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        arg=getNextArgument(arg,true);
        oldVolume=volume;
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        arg5=0;
        loadedConfiguration=0;
    }
    printf("\nTime for equilibrations: %ldsec ,  time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n");
}
