#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "MTGenerator.h"

/*State of optional functions:  (remember to always update text below after changing this file)
 no programmed functions
*/

int N,N2,gaps,activeN,activeN2,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,ODFLength,OCFMode,skipFirstIteration,saveConfigurations,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    onlyMath[2]={0,0},initMode,autoEqLength,neighUpdatingFrequency;
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,intervalOrientations,savedConfigurationsInt,generatorStartPoint=0;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg[2],deltaR,deltaPhi,deltaV=0.1,deltaVND=0.1,
       TReduced,invPotPower,randomStartStep[2],
       neighRadius,neighRadius2,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][3],pi=M_PI,polydisperseDelta,
       boxMatrixInnerParameters[11]; //innerParameters: detBoxMatrix,volume,bm1221m1122,bm0122m0221,bm0211m0112,bm1220m1022,bm0022m0220,bm0210m0012,bm1120m1021,bm0021m0120,bm0110m0011;
double VcpPerParticle,dr[4];
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferFolderIndex[5],bufferTReduced[20],bufferPolDisDel[20],bufferInvPotPow[20],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     orientationsFileName[200]="Orientations",
     orientatCorrelFunFileName[200]="OrientatCorrelFun",
     orientationsResultsFileName[200]="OrientatRes",
     loadConfigurationsFileName[200]="Configurations",
     loadedJOBID[50]="j-none";
/////////////////  PARTICLE functions{
typedef struct {
    double r[3], normR[3];  //globalnie potrzebne przy inicjalizacji, a następnie pełnią rolę współrzędnych względnych
    int bond, index,  //potrzebne do inicjalizacji, a następnie, po wygenerowaniu dimerów zmienna 'bond' wskazuje na indeks dimera w tablicy globalnej dimerów
        neighbours[50], neighCounter;  //uwaga na przepełnienie neighbours[] - wówczas przy tworzeniu listy sąsiadów program może 'pisać' w innych zmiennych, np. w r[] zamiast informować (już tak bywało)
    double diameter, neighEnergyFactor[50];  //index dysku w tablicy globalnej dysków
} atom;

typedef struct {
    double r[3], normR[3], phi[3]; //współrzędne i wersor orientacji wzdłuż osi dimera, potrzebne do przesuwania i obracania 'całego' dimera
    double halfSigma;       //półoś dimeru, od środka do jednego dysku
    atom *atoms[2];         //wskaźniki do obiektów dysków, zastępują ich pole 'bond' wiążąc je w dimer
} particle;


void updateBoxMatrixParameters (double box[3][3], double boxInnerParameters[11]) {
    boxInnerParameters[0]=box[0][0]*box[1][1]*box[2][2]+box[0][1]*box[1][2]*box[2][0]+
                          box[0][2]*box[1][0]*box[2][1]-box[0][2]*box[1][1]*box[2][0]-
                          box[0][0]*box[1][2]*box[2][1]-box[0][1]*box[1][0]*box[2][2];
    boxInnerParameters[1]=fabs(boxInnerParameters[0]);
    boxInnerParameters[2]=box[1][2]*box[2][1]-box[1][1]*box[2][2];
    boxInnerParameters[3]=box[0][1]*box[2][2]-box[0][2]*box[2][1];
    boxInnerParameters[4]=box[0][2]*box[1][1]-box[0][1]*box[1][2];
    boxInnerParameters[5]=box[1][2]*box[2][0]-box[1][0]*box[2][2];
    boxInnerParameters[6]=box[0][0]*box[2][2]-box[0][2]*box[2][0];
    boxInnerParameters[7]=box[0][2]*box[1][0]-box[0][0]*box[1][2];
    boxInnerParameters[8]=box[1][1]*box[2][0]-box[1][0]*box[2][1];
    boxInnerParameters[9]=box[0][0]*box[2][1]-box[0][1]*box[2][0];
    boxInnerParameters[10]=box[0][1]*box[1][0]-box[0][0]*box[1][1];
}

int getIndexFromList (int *list, int listLength, int searchedElement) {
    int indexInList=-1;
    for (int i=0;i<listLength;i++) if (list[i]==searchedElement) {
        indexInList=i;
        break;
    }
    return indexInList;
}

void getAtomsDistanceSquared (atom *a1, atom *a2, particle *p1, particle *p2, double boxMatrix[3][3], bool relativeDiscsInputCoordinates) {
    double normalizedDR[3],roundedNormalizedDR[3],DR[3];
    for (int i=0;i<3;i++) {
        normalizedDR[i]=a1->normR[i]-a2->normR[i];
        DR[i]=a1->r[i]-a2->r[i];
    }
    if (relativeDiscsInputCoordinates) for (int i=0;i<3;i++) {
        normalizedDR[i]+=p1->normR[i]-p2->normR[i];
        DR[i]+=p1->r[i]-p2->r[i];
    }
    for (int i=0;i<3;i++) roundedNormalizedDR[i]=round(normalizedDR[i]);
    for (int j=0;j<3;j++) for (int i=0;i<3;i++) DR[i]-=roundedNormalizedDR[j]*boxMatrix[i][j];
    dr[3]=0; for (int i=0;i<3;i++) {dr[i]=DR[i]; dr[3]+=dr[i]*dr[i];}
}

void adjustNeighRadius (double volume) {
    neighRadius=1.4*cbrt(volume)/cbrt(N2);
    neighRadius2=neighRadius*neighRadius;
}

void updateAtomsNeighbourList (atom *atoms, particle *particles, double boxMatrix[3][3], double volume, bool relativeDiscsInputCoordinates) {
    adjustNeighRadius(volume);
    for (int i=0;i<activeN2;i++) atoms[i].neighCounter=0;
    for (int i=0;i<activeN2-1;i++) for (int j=i+1;j<activeN2;j++) {
        getAtomsDistanceSquared(&atoms[i],&atoms[j],&particles[atoms[i].bond],&particles[atoms[j].bond],boxMatrix,relativeDiscsInputCoordinates);
        if (dr[3]<neighRadius2 && (atoms[i].bond!=atoms[j].bond || !relativeDiscsInputCoordinates)) {
            atoms[i].neighbours[atoms[i].neighCounter++]=j;
            atoms[j].neighbours[atoms[j].neighCounter++]=i;
        }
    }
}

void computeAtomsPositions (particle *p, double boxParameters[11]) {
    for (int i=0;i<3;i++) {
        p->atoms[0]->r[i]=p->halfSigma*p->phi[i];
        p->atoms[1]->r[i]=-p->halfSigma*p->phi[i];
    }
    for (int i=0;i<2;i++) {
        p->atoms[i]->normR[0]=-(p->atoms[i]->r[0]*boxParameters[2]+p->atoms[i]->r[1]*boxParameters[3]+p->atoms[i]->r[2]*boxParameters[4])/boxParameters[0];
        p->atoms[i]->normR[1]=(p->atoms[i]->r[0]*boxParameters[5]+p->atoms[i]->r[1]*boxParameters[6]+p->atoms[i]->r[2]*boxParameters[7])/boxParameters[0];
        p->atoms[i]->normR[2]=-(p->atoms[i]->r[0]*boxParameters[8]+p->atoms[i]->r[1]*boxParameters[9]+p->atoms[i]->r[2]*boxParameters[10])/boxParameters[0];
    }
}

void checkSinglePeriodicBoundaryConditions (particle *p, double boxMatrix[3][3]) {
    for (int j=0;j<3;j++) {
        for (int i=0;i<3;i++) p->r[i]-=floor(p->normR[j])*boxMatrix[i][j];
        p->normR[j]=fmod(p->normR[j],1); if (p->normR[j]<0) p->normR[j]++;
    }
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[3][3]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

int createRandomGaps (particle *particles, double boxMatrix[3][3], double volume) { //TODO: 1) głupie swapowanie - po co przesuwać cząstki w liście, jak wystarczy po prostu 'przerzucać' gapy na (aktualny) koniec listy i zmniejszać jej wielkość (ucinając automatycznie ten ostatni element). Pamiętać o updatowaniu indeksu dimeru w obiektach dysków.; 2) nie zamienia wszystkich pól obiektów cząstek (przestarzałe); 3) należy również przerzucać na tył listy obiekty dysków przyłączonych do dimerów
    printf("Creating %d random gaps... ",gaps); fflush(stdout);
    adjustNeighRadius(volume);
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTRandom0to1(randomStartStep)*N);
            for (int j=0;j<i;j++) {
                for (int k=0;k<2;k++) for (int l=0;l<2;l++) {
                    getAtomsDistanceSquared(particles[gapsIndexes[i]].atoms[k],particles[gapsIndexes[j]].atoms[l],&particles[gapsIndexes[i]],&particles[gapsIndexes[j]],boxMatrix,true);
                    if (dr[3]<neighRadius2) {
                        i--;
                        innerAttempt++;
                        break;
                    }
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
            for (int j=0;j<3;j++) {
                particles[i].r[j]=particles[i+actualGapIndex].r[j];
                particles[i].normR[j]=particles[i+actualGapIndex].normR[j];
                //TODO: brakujace pola do przepisania
            }
        }
        printf("done\n"); fflush(stdout);
        return 1;
    }
}

double getInitDiameterByDistributionFunction (double (*function)(double), double randomMin, double randomMax) {  //funkcja(double) MUSI być przeskalowana tak, by jej (gęstości prawdopobieństwa) MAX był równy 1
    double diameter;
    do diameter=randomMin+(randomMax-randomMin)*MTRandom0to1(randomStartStep); while (MTRandom0to1(randomStartStep)>function(diameter));
    return diameter;
}

double gaussianFunction (double x) {
    double my=1, delta=polydisperseDelta;  //my-wartosc srednia rozkladu, delta-odchylenie standardowe
    return /* 1/delta/sqrt(2*M_PI)* */exp(-pow(x-my,2)/(2*delta*delta));  //wykomentowany człon to amplituda funkcji w punkcie MAX - w taki sposób funkcja jest skalowana do maxValue=1
}

int initPositions (atom *atoms, particle *particles, double boxMatrix[3][3], double boxParameters[11], double matrixOfParticlesSize[3], int n[3], double matrixCellXYZ[6][6][6][3], double pacFrac) {
    if (generatorStartPoint==0) {
        printf("Setting start position of p-random number generator to actual CPU time (for INIT PROCEDURES)...\n");
        InitRandomMT();
    } else {
        printf("Setting start position of p-random number generator to %ld (for INIT PROCEDURES)...\n",generatorStartPoint);
        InitMT((unsigned int)generatorStartPoint);
    }

    //stage #1 - ustawienie atomów w sieci fcc (face-centered cubic)
    double mod=cbrt(N2/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]/matrixOfParticlesSize[2]), interval[3][3], actualPosition[3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) interval[i][j]=boxMatrix[i][j]/matrixOfParticlesSize[j]/mod*n[j];
    int layerCounter=0, rowCounter=0, columnCounter=0;
    for (int i=0;i<N2;i++) {
        int cellNumber[3][3]={{columnCounter/n[0],rowCounter/n[1],layerCounter/n[2]},{columnCounter%n[0],rowCounter%n[1],layerCounter%n[2]}}; //cellNumber[0/1][X/Y/Z]: 0-numer komorki, 1-kolumna/rzad W komorce
        for (int j=0;j<3;j++) actualPosition[j]=cellNumber[0][0]*interval[j][0]+cellNumber[0][1]*interval[j][1]+cellNumber[0][2]*interval[j][2]+matrixCellXYZ[cellNumber[1][0]][cellNumber[1][1]][cellNumber[1][2]][j]*cbrt(pacFrac);

        for (int j=0;j<3;j++) atoms[i].r[j]=actualPosition[j];
        atoms[i].normR[0]=-(atoms[i].r[0]*boxParameters[2]+atoms[i].r[1]*boxParameters[3]+atoms[i].r[2]*boxParameters[4])/boxParameters[0];
        atoms[i].normR[1]=(atoms[i].r[0]*boxParameters[5]+atoms[i].r[1]*boxParameters[6]+atoms[i].r[2]*boxParameters[7])/boxParameters[0];
        atoms[i].normR[2]=-(atoms[i].r[0]*boxParameters[8]+atoms[i].r[1]*boxParameters[9]+atoms[i].r[2]*boxParameters[10])/boxParameters[0];
        atoms[i].bond=-1; atoms[i].index=-1;

        columnCounter++;
        if (columnCounter*1.000001>=matrixOfParticlesSize[0]*mod) {
            columnCounter=0;
            rowCounter++;
            if (rowCounter*1.000001>=matrixOfParticlesSize[1]*mod) {
                rowCounter=0;
                layerCounter++;
            }
        }
    }

    //stage #2 - tworzenie par sąsiednich dysków (DC crystal  (see NarWojKow2013JCP))
    /*Metoda polega na wylosowaniu 2 losowych NIEzbondowanych cząstek, a następnie skupia się na dwuetapowym postępowaniu dla pierwszej z nich (nazwana: parent).
      Etap (1) - poszukiwanie (losowe) NIEzbondowanej cząstki wśród sąsiadów parenta, a jak się uda, zbondowanie z nią. W przypadku powodzenia, losowane są nowe cząstki parent/target.
      Etap (2) - jeżeli wśród sąsiadów parenta NIE ma NIEzbondowanych cząstek, to wybierany jest z tych sąsiadów (losowo) taki, którego cząstka zbondowana jest bliżej targetu niż parent i
                 wykonuje się 'swapowanie', czyli parent i jego wybrany sąsiad tworzą nowy bond, a 'stary' bond wybranego sąsiada staje się nowym parentem. W efekcie, wylosowane NIEzbondowane
                 cząstki (parent/target) zbliżają się do siebie. W początkowych etapach zbliżanie nie będzie następowało, bo etap (1) będzie prawie zawsze kończył operację, ale gdy bondów będzie
                 już bardzo dużo, to etap (2) będzie coraz częściej wykonywany. Jednocześnie, powinien być on zawsze 'sensownie zbieżny', właśnie przez brak pełnej losowości, a przez losowe, ale jednak
                 'zbliżanie' parentu/targetu (gdy zostaną już tylko 2 NIEzbondowane cząstki, operacja się zakończy, gdy parent/target się zbliżą i połączą bondem).
      Uwaga: w wersji 3D algorytm potrafił wpadać w nieskończoną pętlę swapu, związaną prawdopodobnie z tym, że ze względów numerycznych miał problem z rozstrzygnięciem czy bond sąsiada czy parent
             jest bliżej targetu. Potrafił wybierać wówczas bonda sąsiada, a potem, rozważając etap (2) dla bonda, wybierał znowu 'poprzedniego parenta' i tak się blokował przeskakując między tymi
             dwiema cząstkami, które w rzeczywistości były po prostu w TEJ SAMEJ odległości, a jednocześnie nie istnieli sąsiedzi z bondami bliższymi targetu (bo losowe wybieranie sąsiadów w takiej
             sytuacji w końcu by sprawę załatwiło). W takiej sytuacji, przy braku bondu i swapu nie generuje się błędu, tylko generuje się zupełnie losowy swap, niekoniecznie z cząstką posiadającą
             bond bliższy targetu. Brak jest również '=' w warunku (...)>=dr[2] przy poszukiwaniu swapu, co prawdopodobnie nie ma już jednak znaczenia.
    */
    updateAtomsNeighbourList(atoms,NULL,boxMatrix,boxParameters[1],false);
    int nonAssignedDiscsCount=N2, nonAssignedDiscsList[nonAssignedDiscsCount];
    for (int i=0;i<nonAssignedDiscsCount;i++) nonAssignedDiscsList[i]=i;
    while (nonAssignedDiscsCount>0) {
        int selectedIndex[2]={(int)(MTRandom0to1(randomStartStep)*nonAssignedDiscsCount),(int)(MTRandom0to1(randomStartStep)*nonAssignedDiscsCount)};
        while (selectedIndex[0]==selectedIndex[1]) selectedIndex[1]=(int)(MTRandom0to1(randomStartStep)*nonAssignedDiscsCount);
        //printf("  C1: pInd:%d, tInd:%d, nADCount:%d\n",nonAssignedDiscsList[selectedIndex[0]],nonAssignedDiscsList[selectedIndex[1]],nonAssignedDiscsCount);
        bool bond=false;
        while (!bond) {
            //próba znalezienia niesparowanego, wśród sąsiadów aktualnie rozważanego 'niesparowanego defektu (selectedIndex[0])'
            int parentIndex=nonAssignedDiscsList[selectedIndex[0]], targetIndex=nonAssignedDiscsList[selectedIndex[1]],
                actIndex, randStart=actIndex=(int)(MTRandom0to1(randomStartStep)*atoms[parentIndex].neighCounter),
                direction=MTRandom0to1(randomStartStep)<0.5?-1:1;
            do {
                if (atoms[atoms[parentIndex].neighbours[actIndex]].bond==-1) {
                    //printf("  C2: bond generated:%d->%d\n",parentIndex,atoms[parentIndex].neighbours[actIndex]);
                    atoms[parentIndex].bond=atoms[parentIndex].neighbours[actIndex];
                    atoms[atoms[parentIndex].neighbours[actIndex]].bond=parentIndex;
                    nonAssignedDiscsList[getIndexFromList(nonAssignedDiscsList,nonAssignedDiscsCount,atoms[parentIndex].neighbours[actIndex])]=nonAssignedDiscsList[nonAssignedDiscsCount-1]; nonAssignedDiscsCount--;
                    nonAssignedDiscsList[getIndexFromList(nonAssignedDiscsList,nonAssignedDiscsCount,parentIndex)]=nonAssignedDiscsList[nonAssignedDiscsCount-1]; nonAssignedDiscsCount--; //wyszukiwanie na wypadek, gdyby parentIndex==lastInList, wówczas jego pozycja (selectedIndex[0]) została zmieniona (na poprzednią pozycję podłączanego sąsiada)
                    bond=true; break;
                }
                actIndex+=direction;
                actIndex=actIndex>=atoms[parentIndex].neighCounter?0:actIndex<0?atoms[parentIndex].neighCounter-1:actIndex;
            } while (actIndex!=randStart);
            //przemieszczenie 'niesparowanego defektu (selectedIndex[0])' w kierunku 'targetu (selectedIndex[1])'
            bool swap=false;
            if (!bond) {
                getAtomsDistanceSquared(&atoms[parentIndex],&atoms[targetIndex],NULL,NULL,boxMatrix,false);
                double drSquaredToTarget=dr[3];
                do {
                    getAtomsDistanceSquared(&atoms[targetIndex],&atoms[atoms[atoms[parentIndex].neighbours[actIndex]].bond],NULL,NULL,boxMatrix,false);
                    //printf("  C3: PtoT(%d->%d):%.3E, AtoT(%d[%d-%d]):%.3E\n",parentIndex,targetIndex,sqrt(drSquaredToTarget),atoms[parentIndex].neighbours[actIndex],actIndex,atoms[atoms[parentIndex].neighbours[actIndex]].bond,sqrt(dr[3]));
                    if (drSquaredToTarget>dr[3]*1.000001) { //mnożenie ze względu na 'numerykę' - istotne, aby swap był dokonywany JEDYNIE gdy cząstka bondowa jest rzeczywiście bliżej, a nie w tej samej odległości (ale numerycznie ciut większej), aby uniknąć nieskończonej pętli wymieniania
                        //printf("  C4a: swap generated:%d-%d->%d-%d, nADCount:%d\n",atoms[parentIndex].neighbours[actIndex],atoms[atoms[parentIndex].neighbours[actIndex]].bond,atoms[parentIndex].neighbours[actIndex],parentIndex,nonAssignedDiscsCount);
                        int buffer=atoms[atoms[parentIndex].neighbours[actIndex]].bond;
                        atoms[parentIndex].bond=atoms[parentIndex].neighbours[actIndex];
                        atoms[atoms[parentIndex].neighbours[actIndex]].bond=parentIndex;
                        atoms[buffer].bond=-1;
                        nonAssignedDiscsList[selectedIndex[0]]=buffer;
                        swap=true; break;
                    }
                    actIndex+=direction;
                    actIndex=actIndex>=atoms[parentIndex].neighCounter?0:actIndex<0?atoms[parentIndex].neighCounter-1:actIndex;
                } while (actIndex!=randStart);
            }
            if (!bond && !swap) {//czasami może nastąpić taka 'blokada' - wówczas dokonuje się po prostu swapa z losowym sąsiadem (nie z tym, który ma najbliższego bonda do targetu, bo taki warunek może wpaść w nieskończoną pętlę wymieniania się między tymi dwoma cząstkami)
                int randSwap=(int)(MTRandom0to1(randomStartStep)*atoms[parentIndex].neighCounter);
                //printf("  C4b: random swap generated:%d-%d->%d-%d, nADCount:%d\n",atoms[parentIndex].neighbours[randSwap],atoms[atoms[parentIndex].neighbours[randSwap]].bond,atoms[parentIndex].neighbours[randSwap],parentIndex,nonAssignedDiscsCount);
                int buffer=atoms[atoms[parentIndex].neighbours[randSwap]].bond;
                atoms[parentIndex].bond=atoms[parentIndex].neighbours[randSwap];
                atoms[atoms[parentIndex].neighbours[randSwap]].bond=parentIndex;
                atoms[buffer].bond=-1;
                nonAssignedDiscsList[selectedIndex[0]]=buffer;
            }
        }
    }

    //stage #3 - generowanie dimerów na podstawie par dysków
    int startDiscIndex=0;
    for (int i=0;i<N;i++) {
        for (int j=startDiscIndex;j<N2;j++) if (atoms[j].index==-1) {
            particles[i].atoms[0]=&atoms[j]; particles[i].atoms[1]=&atoms[atoms[j].bond];
            particles[i].atoms[0]->index=j; particles[i].atoms[1]->index=atoms[j].bond;
            for (int k=0;k<2;k++) particles[i].atoms[k]->bond=i;
            startDiscIndex=j; break;
        }
        double roundedNormalizedDR[3];
        for (int j=0;j<3;j++) roundedNormalizedDR[j]=round(particles[i].atoms[0]->normR[j]-particles[i].atoms[1]->normR[j]);
        for (int j=0;j<3;j++) particles[i].r[j]=(particles[i].atoms[0]->r[j]-(roundedNormalizedDR[0]*boxMatrix[j][0]+roundedNormalizedDR[1]*boxMatrix[j][1]+roundedNormalizedDR[2]*boxMatrix[j][2])+particles[i].atoms[1]->r[j])*0.5;
        particles[i].normR[0]=-(particles[i].r[0]*boxParameters[2]+particles[i].r[1]*boxParameters[3]+particles[i].r[2]*boxParameters[4])/boxParameters[0];
        particles[i].normR[1]=(particles[i].r[0]*boxParameters[5]+particles[i].r[1]*boxParameters[6]+particles[i].r[2]*boxParameters[7])/boxParameters[0];
        particles[i].normR[2]=-(particles[i].r[0]*boxParameters[8]+particles[i].r[1]*boxParameters[9]+particles[i].r[2]*boxParameters[10])/boxParameters[0];
        for (int j=0;j<3;j++) particles[i].phi[j]=(particles[i].atoms[0]->r[j]-particles[i].atoms[1]->r[j]-(roundedNormalizedDR[0]*boxMatrix[j][0]+roundedNormalizedDR[1]*boxMatrix[j][1]+roundedNormalizedDR[2]*boxMatrix[j][2]))/cbrt(pacFrac);

    }
    checkPeriodicBoundaryConditions(particles,boxMatrix);

    //stage #4 - generowanie polidyspersji
    for (int i=0;i<N;i++) {
        particles[i].halfSigma=0.5;
        switch (initMode) {
            case 0:{    //brak polidyspersji
                for (int j=0;j<2;j++)
                    particles[i].atoms[j]->diameter=1;
            } break;
            case 1:{    //struktura P1 (see NarWojKow2013JCP)
                for (int j=0;j<2;j++)
                    particles[i].atoms[j]->diameter=getInitDiameterByDistributionFunction(gaussianFunction,1-6*polydisperseDelta,1+6*polydisperseDelta);
            } break;
            case 2:{    //struktura P2 (see NarWojKow2013JCP)
                particles[i].atoms[0]->diameter=getInitDiameterByDistributionFunction(gaussianFunction,1-6*polydisperseDelta,1+6*polydisperseDelta);
                particles[i].atoms[1]->diameter=2-particles[i].atoms[0]->diameter;
            } break;
            case 3:{    //struktura P3 (see NarWojKow2013JCP)
                for (int j=0;j<2;j++)
                    particles[i].atoms[j]->diameter=getInitDiameterByDistributionFunction(gaussianFunction,1-6*polydisperseDelta,1+6*polydisperseDelta);
                particles[i].halfSigma=(particles[i].atoms[0]->diameter+particles[i].atoms[1]->diameter)*0.25;
            } break;
            case 4:{    //struktura B1 (see NarWojKow2013JCP)
                for (int j=0;j<2;j++)
                    if (MTRandom0to1(randomStartStep)<0.5) particles[i].atoms[j]->diameter=1-polydisperseDelta;
                    else particles[i].atoms[j]->diameter=1+polydisperseDelta;
            } break;
            case 5:{    //struktura B2 (see NarWojKow2013JCP)
                if (MTRandom0to1(randomStartStep)<0.5) particles[i].atoms[0]->diameter=1-polydisperseDelta; else particles[i].atoms[0]->diameter=1+polydisperseDelta;
                particles[i].atoms[1]->diameter=2-particles[i].atoms[0]->diameter;
            } break;
        }
        if (initMode==4) {  //doprowadzenie układu binarnego small/large do udziału po 50% każdego typu (potrzebne tylko dla struktury B1)
            int smallCounter=0;
            for (int i=0;i<N2;i++) if (atoms[i].diameter<1) smallCounter++;
            while (smallCounter>N) {
                int randIndex = (int)(MTRandom0to1(randomStartStep)*N2);
                if (atoms[randIndex].diameter<1) {atoms[randIndex].diameter=1+polydisperseDelta; smallCounter--;}
            }
            while (smallCounter<N) {
                int randIndex = (int)(MTRandom0to1(randomStartStep)*N2);
                if (atoms[randIndex].diameter>1) {atoms[randIndex].diameter=1-polydisperseDelta; smallCounter++;}
            }
        }
        computeAtomsPositions(&particles[i],boxMatrixInnerParameters);
    }
    if (gaps>0) {
        //updateAtomsNeighbourList(atoms,particles,boxMatrix,boxParameters[1],true);  //na ten moment NIE wiem z jakiego powodu kiedyś stwierdziłem, że to jest potrzebne do generowania luk (wydaje mi się, że jednak nie jest...)
        return createRandomGaps(particles,boxMatrix,boxParameters[1]);
    } else return 1;
}

double getEnergyAll (atom *atoms, particle *particles, double boxMatrix[3][3]) {
    double energy=0;
    for (int i=0;i<activeN2;i++) for (int j=0;j<atoms[i].neighCounter;j++) {
        if (i<atoms[i].neighbours[j]) {
            getAtomsDistanceSquared(&atoms[i],&atoms[atoms[i].neighbours[j]],&particles[atoms[i].bond],&particles[atoms[atoms[i].neighbours[j]].bond],boxMatrix,true);
            atoms[i].neighEnergyFactor[j]=pow(0.5*(atoms[i].diameter+atoms[atoms[i].neighbours[j]].diameter)/sqrt(dr[3]),invPotPower);
            atoms[atoms[i].neighbours[j]].neighEnergyFactor[getIndexFromList(atoms[atoms[i].neighbours[j]].neighbours,atoms[atoms[i].neighbours[j]].neighCounter,i)]=atoms[i].neighEnergyFactor[j];
            energy+=atoms[i].neighEnergyFactor[j];
        }
    }
    return energy;
}

double getVicinityEnergyChange (atom *atoms, particle *particles, particle *dispPart, particle *parentPart, double boxMatrix[3][3]) {
    double deltaEnergy=0;
    for (int i=0;i<2;i++) for (int j=0;j<parentPart->atoms[i]->neighCounter;j++) {
        getAtomsDistanceSquared(&atoms[parentPart->atoms[i]->neighbours[j]],dispPart->atoms[i],&particles[atoms[parentPart->atoms[i]->neighbours[j]].bond],dispPart,boxMatrix,true);
        dispPart->atoms[i]->neighEnergyFactor[j]=pow(0.5*(atoms[parentPart->atoms[i]->neighbours[j]].diameter+parentPart->atoms[i]->diameter)/sqrt(dr[3]),invPotPower);
        deltaEnergy+=dispPart->atoms[i]->neighEnergyFactor[j]-parentPart->atoms[i]->neighEnergyFactor[j];
    }
    return deltaEnergy;
}

void generateRandomUnitVector (double randVector[3]) {
    double randSquared,rand1,rand2;
    do {
        rand1=1-2*MTRandom0to1(randomStartStep);
        rand2=1-2*MTRandom0to1(randomStartStep);
        randSquared=rand1*rand1+rand2*rand2;
    } while (randSquared>1);
    double randH=2*sqrt(1-randSquared);
    randVector[0]=rand1*randH; randVector[1]=rand2*randH; randVector[2]=1-2*randSquared;
}

int attemptToDisplaceAParticle (atom *atoms, particle *particles, int index, double boxMatrix[3][3], double boxParameters[11], double *totalInteractionEnergy, bool moveType) {
    int result=1;
    particle displacedParticle;
    displacedParticle.halfSigma=particles[index].halfSigma;
    atom dispAtom[2]; for (int i=0;i<2;i++) displacedParticle.atoms[i]=&dispAtom[i];
    for (int i=0;i<3;i++) {
        displacedParticle.r[i]=particles[index].r[i];
        displacedParticle.phi[i]=particles[index].phi[i];
    }
    if (moveType) { //translational move
        for (int i=0;i<3;i++) displacedParticle.r[i]+=(MTRandom0to1(randomStartStep)-0.5)*deltaR;
        displacedParticle.normR[0]=-(displacedParticle.r[0]*boxParameters[2]+displacedParticle.r[1]*boxParameters[3]+displacedParticle.r[2]*boxParameters[4])/boxParameters[0];
        displacedParticle.normR[1]=(displacedParticle.r[0]*boxParameters[5]+displacedParticle.r[1]*boxParameters[6]+displacedParticle.r[2]*boxParameters[7])/boxParameters[0];
        displacedParticle.normR[2]=-(displacedParticle.r[0]*boxParameters[8]+displacedParticle.r[1]*boxParameters[9]+displacedParticle.r[2]*boxParameters[10])/boxParameters[0];
        for (int i=0;i<2;i++) for (int j=0;j<3;j++) {
            displacedParticle.atoms[i]->r[j]=particles[index].atoms[i]->r[j];
            displacedParticle.atoms[i]->normR[j]=particles[index].atoms[i]->normR[j];
        }
    } else { //orientational move
        for (int i=0;i<3;i++) displacedParticle.normR[i]=particles[index].normR[i];
        double randUnitVector[3],phiLengthInv=0; generateRandomUnitVector(randUnitVector);
        for (int i=0;i<3;i++) {
            displacedParticle.phi[i]+=randUnitVector[i]*deltaPhi;
            phiLengthInv+=displacedParticle.phi[i]*displacedParticle.phi[i];
        } phiLengthInv=1/sqrt(phiLengthInv); for (int i=0;i<3;i++) displacedParticle.phi[i]*=phiLengthInv;
        computeAtomsPositions(&displacedParticle,boxParameters);
    }
    double deltaEnergy=getVicinityEnergyChange(atoms,particles,&displacedParticle,&particles[index],boxMatrix);

    double arg=-deltaEnergy/TReduced; //T*=kT/eps [eps - jednostka energii (w domyśle =1)]
    if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
    else {
        *totalInteractionEnergy+=deltaEnergy;
        if (moveType) {
            for (int i=0;i<3;i++) {
                particles[index].r[i]=displacedParticle.r[i];
                particles[index].normR[i]=displacedParticle.normR[i];
            }
            checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix);
        } else for (int i=0;i<3;i++) {
            particles[index].phi[i]=displacedParticle.phi[i];
            for (int j=0;j<2;j++) {
                particles[index].atoms[j]->r[i]=displacedParticle.atoms[j]->r[i];
                particles[index].atoms[j]->normR[i]=displacedParticle.atoms[j]->normR[i];
            }
        }
        for (int i=0;i<2;i++) for (int j=0;j<particles[index].atoms[i]->neighCounter;j++) {
            particles[index].atoms[i]->neighEnergyFactor[j]=displacedParticle.atoms[i]->neighEnergyFactor[j];
            atoms[particles[index].atoms[i]->neighbours[j]].neighEnergyFactor[getIndexFromList(atoms[particles[index].atoms[i]->neighbours[j]].neighbours,atoms[particles[index].atoms[i]->neighbours[j]].neighCounter,particles[index].atoms[i]->index)]=displacedParticle.atoms[i]->neighEnergyFactor[j];
        }
    }
    return result;
}

void cloneParticlesForSpecificBoxMatrix (particle *clonedParticles, particle *particles, atom *clonedAtoms, atom *atoms, double boxMatrix[3][3], double boxParameters[11]) {
    for (int i=0;i<activeN;i++) {
        for (int j=0;j<3;j++) clonedParticles[i].normR[j]=particles[i].normR[j];
        for (int j=0;j<3;j++) clonedParticles[i].r[j]=boxMatrix[j][0]*particles[i].normR[0]+boxMatrix[j][1]*particles[i].normR[1]+boxMatrix[j][2]*particles[i].normR[2];
        for (int j=0;j<2;j++) clonedParticles[i].atoms[j]=&clonedAtoms[particles[i].atoms[j]->index];
    }
    for (int i=0;i<activeN2;i++) {
        for (int j=0;j<3;j++) clonedAtoms[i].r[j]=atoms[i].r[j];
        clonedAtoms[i].normR[0]=-(atoms[i].r[0]*boxParameters[2]+atoms[i].r[1]*boxParameters[3]+atoms[i].r[2]*boxParameters[4])/boxParameters[0];
        clonedAtoms[i].normR[1]=(atoms[i].r[0]*boxParameters[5]+atoms[i].r[1]*boxParameters[6]+atoms[i].r[2]*boxParameters[7])/boxParameters[0];
        clonedAtoms[i].normR[2]=-(atoms[i].r[0]*boxParameters[8]+atoms[i].r[1]*boxParameters[9]+atoms[i].r[2]*boxParameters[10])/boxParameters[0];
    }
}

int attemptToChangeVolume (atom *atoms, particle *particles, double pressure, double boxMatrix[3][3], double boxParameters[11], double *totalInteractionEnergy) {
    int result=1;
    double newBoxMatrix[3][3], lnNewBoxMatrix[3], newBoxMatrixInnerParameters[11], newTotalInteractionEnergy=0;
    if (boxParameters[1]/VcpPerParticle/N2<pressureRealOfNotFluid) {
    //if (pressure>pressureRealOfNotFluid) {  //dozwolone zmiany ksztaltu pudla (faza stala)
        //ruchy liniowe dla wszystkich elementów pudła
        /*newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[1][1]=boxMatrix[1][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[2][2]=boxMatrix[2][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][2]=boxMatrix[0][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[1][2]=boxMatrix[1][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;*/

        //ruchy logarytmiczne dla elementów diagonalnych oraz liniowe dla niediagonalnych (logarytmiczne lepiej próbkują przestrzeń dla stanu wyjściowego >1, a liniowe lepiej dla stanów <1 (szczególnie w okolicy 0))
        for (int i=0;i<3;i++) lnNewBoxMatrix[i]=log(boxMatrix[i][i])+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        for (int i=0;i<3;i++) newBoxMatrix[i][i]=exp(lnNewBoxMatrix[i]);
        newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][2]=boxMatrix[0][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[1][2]=boxMatrix[1][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        double modifier=newBoxMatrix[0][0]/boxMatrix[0][0];
        newBoxMatrix[1][1]=boxMatrix[1][1]*modifier;
        newBoxMatrix[2][2]=boxMatrix[2][2]*modifier;
        newBoxMatrix[0][1]=boxMatrix[0][1]*modifier;
        newBoxMatrix[0][2]=boxMatrix[0][2]*modifier;
        newBoxMatrix[1][2]=boxMatrix[1][2]*modifier;
    }
    newBoxMatrix[1][0]=newBoxMatrix[0][1]; newBoxMatrix[2][0]=newBoxMatrix[0][2]; newBoxMatrix[2][1]=newBoxMatrix[1][2];
    updateBoxMatrixParameters(newBoxMatrix,newBoxMatrixInnerParameters);

    particle particlesInNewBox[activeN]; atom atomsInNewBox[activeN2];
    cloneParticlesForSpecificBoxMatrix(particlesInNewBox,particles,atomsInNewBox,atoms,newBoxMatrix,newBoxMatrixInnerParameters);
    for (int i=0;i<activeN2;i++) for (int j=0;j<atoms[i].neighCounter;j++)
        if (i<atoms[i].neighbours[j]) {
            getAtomsDistanceSquared(&atomsInNewBox[i],&atomsInNewBox[atoms[i].neighbours[j]],&particlesInNewBox[atoms[i].bond],&particlesInNewBox[atoms[atoms[i].neighbours[j]].bond],newBoxMatrix,true);
            atomsInNewBox[i].neighEnergyFactor[j]=pow(0.5*(atoms[i].diameter+atoms[atoms[i].neighbours[j]].diameter)/sqrt(dr[3]),invPotPower);
            newTotalInteractionEnergy+=atomsInNewBox[i].neighEnergyFactor[j];
        }    
    double arg=-((newTotalInteractionEnergy-(*totalInteractionEnergy)+pressure*(newBoxMatrixInnerParameters[1]-boxParameters[1]))/TReduced-(((double)N+1.0)*log(newBoxMatrixInnerParameters[1]/boxParameters[1])+log((newBoxMatrix[0][0]+newBoxMatrix[1][1]+newBoxMatrix[2][2])/(boxMatrix[0][0]+boxMatrix[1][1]+boxMatrix[2][2]))));
    if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
    else {
        for (int i=0;i<11;i++) boxParameters[i]=newBoxMatrixInnerParameters[i];
        *totalInteractionEnergy=newTotalInteractionEnergy;
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
        for (int i=0;i<activeN;i++) for (int j=0;j<3;j++) particles[i].r[j]=particlesInNewBox[i].r[j];
        for (int i=0;i<activeN2;i++) {
            for (int j=0;j<3;j++) atoms[i].normR[j]=atomsInNewBox[i].normR[j];
            for (int j=0;j<atoms[i].neighCounter;j++) if (i<atoms[i].neighbours[j]) {
                atoms[i].neighEnergyFactor[j]=atomsInNewBox[i].neighEnergyFactor[j];
                atoms[atoms[i].neighbours[j]].neighEnergyFactor[getIndexFromList(atoms[atoms[i].neighbours[j]].neighbours,atoms[atoms[i].neighbours[j]].neighCounter,i)]=atomsInNewBox[i].neighEnergyFactor[j];
            }
        }
    }
    return result;
}

int attemptToChangeVolumeSeparate (atom *atoms, particle *particles, double pressure, double boxMatrix[3][3], double boxParameters[11], double *totalInteractionEnergy, double deltaV, bool moveType) {
    int result=1;
    double newBoxMatrix[3][3], newBoxMatrixInnerParameters[11], newTotalInteractionEnergy=0, NFactor;
    if (boxParameters[1]/VcpPerParticle/N2<pressureRealOfNotFluid) {
        if (moveType) { //logarytmiczny ruch elementów diagonalnych
            NFactor=N+1;
            for (int i=0;i<3;i++) newBoxMatrix[i][i]=exp(log(boxMatrix[i][i])+(MTRandom0to1(randomStartStep)-0.5)*deltaV);
            newBoxMatrix[0][1]=boxMatrix[0][1];
            newBoxMatrix[0][2]=boxMatrix[0][2];
            newBoxMatrix[1][2]=boxMatrix[1][2];
        } else {        //liniowy ruch elementów niediagonalnych
            NFactor=N;
            for (int i=0;i<3;i++) newBoxMatrix[i][i]=boxMatrix[i][i];
            newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
            newBoxMatrix[0][2]=boxMatrix[0][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
            newBoxMatrix[1][2]=boxMatrix[1][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        }
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        NFactor=N+1;
        newBoxMatrix[0][0]=exp(log(boxMatrix[0][0])+(MTRandom0to1(randomStartStep)-0.5)*deltaV);
        double modifier=newBoxMatrix[0][0]/boxMatrix[0][0];
        newBoxMatrix[1][1]=boxMatrix[1][1]*modifier;
        newBoxMatrix[2][2]=boxMatrix[2][2]*modifier;
        newBoxMatrix[0][1]=boxMatrix[0][1]*modifier;
        newBoxMatrix[0][2]=boxMatrix[0][2]*modifier;
        newBoxMatrix[1][2]=boxMatrix[1][2]*modifier;
    }
    newBoxMatrix[1][0]=newBoxMatrix[0][1]; newBoxMatrix[2][0]=newBoxMatrix[0][2]; newBoxMatrix[2][1]=newBoxMatrix[1][2];
    updateBoxMatrixParameters(newBoxMatrix,newBoxMatrixInnerParameters);

    particle particlesInNewBox[activeN]; atom atomsInNewBox[activeN2];
    cloneParticlesForSpecificBoxMatrix(particlesInNewBox,particles,atomsInNewBox,atoms,newBoxMatrix,newBoxMatrixInnerParameters);
    for (int i=0;i<activeN2;i++) for (int j=0;j<atoms[i].neighCounter;j++)
        if (i<atoms[i].neighbours[j]) {
            getAtomsDistanceSquared(&atomsInNewBox[i],&atomsInNewBox[atoms[i].neighbours[j]],&particlesInNewBox[atoms[i].bond],&particlesInNewBox[atoms[atoms[i].neighbours[j]].bond],newBoxMatrix,true);
            atomsInNewBox[i].neighEnergyFactor[j]=pow(0.5*(atoms[i].diameter+atoms[atoms[i].neighbours[j]].diameter)/sqrt(dr[3]),invPotPower);
            newTotalInteractionEnergy+=atomsInNewBox[i].neighEnergyFactor[j];
        }
    double arg=-((newTotalInteractionEnergy-(*totalInteractionEnergy)+pressure*(newBoxMatrixInnerParameters[1]-boxParameters[1]))/TReduced-(NFactor*log(newBoxMatrixInnerParameters[1]/boxParameters[1])+log((newBoxMatrix[0][0]+newBoxMatrix[1][1]+newBoxMatrix[2][2])/(boxMatrix[0][0]+boxMatrix[1][1]+boxMatrix[2][2]))));
    if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
    else {
        for (int i=0;i<11;i++) boxParameters[i]=newBoxMatrixInnerParameters[i];
        *totalInteractionEnergy=newTotalInteractionEnergy;
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
        for (int i=0;i<activeN;i++) for (int j=0;j<3;j++) particles[i].r[j]=particlesInNewBox[i].r[j];
        for (int i=0;i<activeN2;i++) {
            for (int j=0;j<3;j++) atoms[i].normR[j]=atomsInNewBox[i].normR[j];
            for (int j=0;j<atoms[i].neighCounter;j++) if (i<atoms[i].neighbours[j]) {
                atoms[i].neighEnergyFactor[j]=atomsInNewBox[i].neighEnergyFactor[j];
                atoms[atoms[i].neighbours[j]].neighEnergyFactor[getIndexFromList(atoms[atoms[i].neighbours[j]].neighbours,atoms[atoms[i].neighbours[j]].neighCounter,i)]=atomsInNewBox[i].neighEnergyFactor[j];
            }
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
    strcpy(buffer,"3D_N-"); strncat(buffer,bufferN,20);
    strncat(buffer,"_gaps-",10); strncat(buffer,bufferGaps,20);
    strncat(buffer,"_G-",5); strncat(buffer,bufferG,5);
    strncat(buffer,"_badanie-",10); strncat(buffer,bufferFolderIndex,5);
    strncat(buffer,"_T-",5); strncat(buffer,bufferTReduced,20);
    strncat(buffer,"_pDD-",6); strncat(buffer,bufferPolDisDel,20);
    strncat(buffer,"_iPP-",6); strncat(buffer,bufferInvPotPow,20);
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

bool isLineCorrect(char linia[4096]) {
    sscanf(linia,"%c",linia);
    int actIndex=0, dataIndex=0; while (dataIndex<6) {
        char data[50]="";
        int licznik=0, dotCounter=0;
        while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
        if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;} //test of single dot in a number
        actIndex++; dataIndex++;
        if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10; //test of dot position after first digit
    } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10; //test of max 6 numbers in a row

    return (dataIndex<10);
}

double getAvErrorFromSumEps (double sum, double denominator) {
    return sqrt(sum/denominator);
}

void updateTableAndGetActualMean (double table[100], double & mean, int const & changeIndex, double const & changeValue) {
    mean-=table[changeIndex]*0.01; table[changeIndex]=changeValue; mean+=changeValue*0.01;
}

int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[500];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,500,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20]="";
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++]; data[licznik]=' ';
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:initMode=strtol(data,NULL,10);break;
                case 4:polydisperseDelta=strtod(data,NULL);break;
                case 5:TReduced=strtod(data,NULL);break;
                case 6:invPotPower=strtod(data,NULL);break;
                case 7:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 8:growing=strtol(data,NULL,10);break;
                case 9:loadedConfiguration=strtol(data,NULL,10);break;
                case 10:loadedArg=strtod(data,NULL);break;
                case 11:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 12:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 13:loadedSetGenerator=strtol(data,NULL,10);break;
                case 14:iterationsNumber=strtol(data,NULL,10);break;
                case 15:intervalSampling=strtol(data,NULL,10);break;
                case 16:intervalOutput=strtol(data,NULL,10);break;
                case 17:saveConfigurations=strtol(data,NULL,10);break;
                case 18:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 19:ODFLength=strtol(data,NULL,10);break;
                case 20:OCFMode=strtol(data,NULL,10);break;
                case 21:neighUpdatingFrequency=strtol(data,NULL,10);break;
                case 22:intervalOrientations=strtol(data,NULL,10);break;
                case 23:skipFirstIteration=strtol(data,NULL,10);break;
                case 24:useSpecificDirectory=strtol(data,NULL,10);break;
                case 25:autoEqLength=strtol(data,NULL,10);break;
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
                    growing=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    TReduced=strtod(arguments[10],NULL);
                    polydisperseDelta=strtod(arguments[11],NULL);
                    invPotPower=strtod(arguments[12],NULL);
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    N=strtol(arguments[3],NULL,10);
                    gaps=strtol(arguments[4],NULL,10);
                    growing=strtol(arguments[5],NULL,10);
                    iterationsNumber=strtol(arguments[6],NULL,10);
                    useSpecificDirectory=strtol(arguments[7],NULL,10);
                    skipFirstIteration=strtol(arguments[8],NULL,10);
                    pointNumber=strtol(arguments[9],NULL,10);
                    generatorStartPoint=strtol(arguments[10],NULL,10);
                    TReduced=strtod(arguments[11],NULL);
                    polydisperseDelta=strtod(arguments[12],NULL);
                    invPotPower=strtod(arguments[13],NULL);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    growing=strtol(arguments[6],NULL,10);
                    loadedConfiguration=1;
                    loadedArg=strtod(arguments[7],NULL);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    TReduced=strtod(arguments[11],NULL);
                    polydisperseDelta=strtod(arguments[12],NULL);
                    invPotPower=strtod(arguments[13],NULL);
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    growing=strtol(arguments[6],NULL,10);
                    loadedConfiguration=1; loadType=1;
                    pointNumber=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    TReduced=strtod(arguments[11],NULL);
                    polydisperseDelta=strtod(arguments[12],NULL);
                    invPotPower=strtod(arguments[13],NULL);
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==13) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    onlyMath[0]=1;
                    onlyMath[1]=strtol(arguments[3],NULL,10);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    growing=strtol(arguments[6],NULL,10);
                    loadedConfiguration=0;
                    pointNumber=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=0;
                    TReduced=strtod(arguments[10],NULL);
                    polydisperseDelta=strtod(arguments[11],NULL);
                    invPotPower=strtod(arguments[12],NULL);
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, growing, iterationsNumber, useSpecificDirectory, TReduced, polydisperseDelta, invPotPower\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, generatorStartPoint, TReduced, polydisperseDelta, invPotPower\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration, TReduced, polydisperseDelta, invPotPower\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, TReduced, polydisperseDelta, invPotPower\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, growing, pointNumber, iterationsNumber, useSpecificDirectory, TReduced, polydisperseDelta, invPotPower\n");
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

    for (int i=0;i<pointNumber;i++) getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg[0];
    N2=N*2; activeN=N-gaps; activeN2=activeN*2;

    if (fabs(cbrt(N2/4.0)-floor(cbrt(N2/4.0)))>0.000001) {
        printf("ERROR: Not supported N: %d.\n",N);
        return 0;
    }

    //stale wynikajace z zadanych parametrow multimerow
    deltaR=deltaPhi=maxDeltaR;  //jednostka odległości to 1

    //nazwy folderow na podstawie parametrow programu
    sprintf(bufferG,"%d",growing); sprintf(bufferN,"%d",N); sprintf(bufferGaps,"%d",gaps); sprintf(bufferTReduced,"%.3E",TReduced);
    sprintf(bufferPolDisDel,"%.3E",polydisperseDelta); sprintf(bufferInvPotPow,"%.3E",invPotPower);

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
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.4E",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(orientationsFileName,JOBID,true);
    if (OCFMode) addAppendix(orientatCorrelFunFileName,JOBID,true);
    addAppendix(orientationsResultsFileName,JOBID,true);

    atom atoms[N2]; particle particles[N];
    double args[15];

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileOrientations, *fileOrientatCorrelFun, *fileAllResults, *fileAllOrientations, *fileOrientationsResults, *fileAllOrientationsResults;
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        fprintf(fileResults,"cycle\tpressure\trho*\tdRho*\tv*\tdV*\tvolume\tdVolume\tH00\tdH00\tH11\tdH11\tH22\tdH22\tH01\tdH01\tH02\tdH02\tH12\tdH12\tB\tdB\tavMy\tdAvMy\tmy1\tdMy1\tmy2\tdMy2\tE\tdE\tnu_100_all\tdNu_100_all\tnu_111_all\tdNu_111_all\tnu_110_1m10\tdNu_110_1m10\tnu_110_001\tdNu_110_001\tS11c\tdS11c\tS12c\tdS12c\tS44c\tdS44c\tS11\tdS11\tS12\tdS12\tS13\tdS13\tS14\tdS14\tS15\tdS15\tS16\tdS16\tS22\tdS22\tS23\tdS23\tS24\tdS24\tS25\tdS25\tS26\tdS26\tS33\tdS33\tS34\tdS34\tS35\tdS35\tS36\tdS36\tS44\tdS44\tS45\tdS45\tS46\tdS46\tS55\tdS55\tS56\tdS56\tS66\tdS66\tB11>0\tB44>0\t-1/2<=B12/B11<=1\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        fprintf(fileExcelResults,"pressure\tv*\tB\tdB\tavMy\tdAvMy\tE\tdE\tnu_100_all\tdNu_100_all\tnu_111_all\tdNu_111_all\tnu_110_1m10\tdNu_110_1m10\tnu_110_001\tdNu_110_001\tS11c\tdS11c\tS12c\tdS12c\tS44c\tdS44c\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg[2]={startArg[0],startArg[1]}, oldBoxMatrix[3][3], oldTotalInteractionEnergy;
    while (growing>=0) {
        double totalInteractionEnergy, pressureReduced=arg[0], pressure=pressureReduced,  //pRed=p*s^3/eps -> p=pRed*eps/s^3, s=1, eps=1 (wymiar p: eps/sigma^3)
               boxMatrix[3][3],matrixOfParticlesSize[3],unitCellAtCP[3],
               matrixCellXYZ[6][6][6][3];
        int n[3]; //n[X/Y/Z], matrixCell[n[XMax]][n[YMax]][n[ZMax]][x/y/z], zatem: n[X/Y/Z](max)=6
        unitCellAtCP[0]=unitCellAtCP[1]=unitCellAtCP[2]=sqrt(2);  //struktura fcc dla kul, jednostka odległości (sigma) to 1
        n[0]=1; n[1]=2; n[2]=2;
        matrixCellXYZ[0][0][0][0]=0; matrixCellXYZ[0][0][0][1]=0; matrixCellXYZ[0][0][0][2]=0;
        matrixCellXYZ[0][1][0][0]=unitCellAtCP[0]/2.0; matrixCellXYZ[0][1][0][1]=unitCellAtCP[1]/2.0; matrixCellXYZ[0][1][0][2]=0;
        matrixCellXYZ[0][0][1][0]=unitCellAtCP[0]/2.0; matrixCellXYZ[0][0][1][1]=0; matrixCellXYZ[0][0][1][2]=unitCellAtCP[2]/2.0;
        matrixCellXYZ[0][1][1][0]=unitCellAtCP[0]; matrixCellXYZ[0][1][1][1]=unitCellAtCP[1]/2.0; matrixCellXYZ[0][1][1][2]=unitCellAtCP[2]/2.0;
        VcpPerParticle=unitCellAtCP[0]*unitCellAtCP[1]*unitCellAtCP[2]/(double)n[0]/(double)n[1]/(double)n[2];

        matrixOfParticlesSize[0]=1; matrixOfParticlesSize[1]=2; matrixOfParticlesSize[2]=2;
        double NLinearMod = cbrt(N2/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]/matrixOfParticlesSize[2]);
        if (startArg[0]==arg[0]) {
            for (int i=0;i<3;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]/(double)n[i]*NLinearMod*(growing?cbrt(startMinPacFrac):cbrt(startMaxPacFrac));
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0; boxMatrix[2][0]=0.0; boxMatrix[0][2]=0.0; boxMatrix[1][2]=0.0; boxMatrix[2][1]=0.0;
        } else {
            for (int i=0;i<3;i++) for (int j=0;j<3;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
            totalInteractionEnergy=oldTotalInteractionEnergy;
        }
        updateBoxMatrixParameters(boxMatrix,boxMatrixInnerParameters);
        double rho=N2/boxMatrixInnerParameters[1], pacFrac=1.0/VcpPerParticle/rho;

        if (!onlyMath[0]) {
            if (arg[0]==startArg[0] && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (StartDen: %.7E, startPacFrac: %.7E, T*: %.3E, invPotPower: %.3E, pDD: %.3E, initMode: %d)\n",N,gaps,growing,startArg[0],rho,pacFrac,TReduced,invPotPower,polydisperseDelta,initMode);
                if (!initPositions(atoms,particles,boxMatrix,boxMatrixInnerParameters,matrixOfParticlesSize,n,matrixCellXYZ,pacFrac)) return 0;
                updateAtomsNeighbourList(atoms,particles,boxMatrix,boxMatrixInnerParameters[1],true);
                for (int i=0;i<N2;i++) if (atoms[i].neighCounter!=11) {
                    printf("ERROR: Wrong number of neighbours (!=11) for at least one (index: %d) atom.\n",i);
                    return 0;
                }
                totalInteractionEnergy=getEnergyAll(atoms,particles,boxMatrix);
            } else if (loadedConfiguration) {
                char configurations[4096];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                for (int i=0;i<3;i++) fgets(configurations,4096,fileCTL);
                int character,dataType=0,pIndex=-1; char data[50]=""; int actIndex=0;
                while (dataType<16) {
                    character=fgetc(fileCTL); //character is in int, but it can work as char
                    if (dataType<15) { //stage #1 configuration parameters
                        if (character==' ') {
                            data[actIndex++]=' '; //end of data without clearing the entire array
                            args[dataType++]=strtod(data,NULL);
                            if (dataType==15) {
                                boxMatrix[0][0]=args[5]; boxMatrix[1][1]=args[6]; boxMatrix[2][2]=args[7];
                                boxMatrix[1][0]=args[8]; boxMatrix[0][1]=args[8];
                                boxMatrix[2][0]=args[9]; boxMatrix[0][2]=args[9];
                                boxMatrix[1][2]=args[10]; boxMatrix[2][1]=args[10];
                                deltaR=args[11]; deltaPhi=args[12]; deltaV=args[13]; deltaVND=args[14];
                                arg[0]=args[3]; pressureReduced=arg[0]; pressure=pressureReduced;

                                updateBoxMatrixParameters(boxMatrix,boxMatrixInnerParameters);
                                rho=N2/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                            }
                            actIndex=0;
                            continue;
                        }
                        data[actIndex++]=character;
                    } else { //stage #2 configuration parameters (coordinates of particles)
                        if (character=='m') {pIndex++; continue;}
                        if (pIndex>=0) {
                            for (int i=0;i<9;i++) {
                                character=fgetc(fileCTL);
                                while (character!=',' && character!=']') {
                                    data[actIndex++]=character;
                                    character=fgetc(fileCTL);
                                } data[actIndex++]=' '; actIndex=0;
                                if (i<3) particles[pIndex].r[i]=strtod(data,NULL);
                                else if (i<6) particles[pIndex].phi[i-3]=strtod(data,NULL);
                                else if (i<7) particles[pIndex].halfSigma=strtod(data,NULL);
                                else {
                                    int discGlobalIndex=pIndex*2+i-7;
                                    particles[pIndex].atoms[i-7]=&atoms[discGlobalIndex];
                                    atoms[discGlobalIndex].bond=pIndex; atoms[discGlobalIndex].index=discGlobalIndex;
                                    atoms[discGlobalIndex].diameter=strtod(data,NULL);
                                }
                            }                            
                            particles[pIndex].normR[0]=-(particles[pIndex].r[0]*boxMatrixInnerParameters[2]+particles[pIndex].r[1]*boxMatrixInnerParameters[3]+particles[pIndex].r[2]*boxMatrixInnerParameters[4])/boxMatrixInnerParameters[0];
                            particles[pIndex].normR[1]=(particles[pIndex].r[0]*boxMatrixInnerParameters[5]+particles[pIndex].r[1]*boxMatrixInnerParameters[6]+particles[pIndex].r[2]*boxMatrixInnerParameters[7])/boxMatrixInnerParameters[0];
                            particles[pIndex].normR[2]=-(particles[pIndex].r[0]*boxMatrixInnerParameters[8]+particles[pIndex].r[1]*boxMatrixInnerParameters[9]+particles[pIndex].r[2]*boxMatrixInnerParameters[10])/boxMatrixInnerParameters[0];

                            computeAtomsPositions(&particles[pIndex],boxMatrixInnerParameters);
                            while (character!=']') character=fgetc(fileCTL); fgetc(fileCTL); //next to read: 'm'
                            if (pIndex>=activeN-1) dataType++;
                        }
                    }
                }
                fclose(fileCTL);
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (startDen: %.7E, startPacFrac: %.7E, T*: %.3E, invPotPower: %.3E, pDD: %.3E), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4E, DeltaPhi: %.4E, DeltaV: %.4E, DeltaVND: %.4E\n",N,gaps,growing,pressureReduced,rho,pacFrac,TReduced,invPotPower,polydisperseDelta,args[0],args[1],(long)args[4],deltaR,deltaPhi,deltaV,deltaVND);
                //for (int i=0;i<3;i++) for (int j=0;j<3;j++) printf("boxMatrix[%d][%d]=%.17E\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<activeN;i++) printf("%d:  %.17E,  %.17E,  %.17E,  %.17E,  %.17E,  %.17E,  %.17E,  %.17E,  %.17E\n",i,particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].phi[0],particles[i].phi[1],particles[i].phi[2],particles[i].halfSigma,particles[i].atoms[0]->diameter,particles[i].atoms[1]->diameter);return 0;
                updateAtomsNeighbourList(atoms,particles,boxMatrix,boxMatrixInnerParameters[1],true);
                totalInteractionEnergy=getEnergyAll(atoms,particles,boxMatrix);
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
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                        if (generatorStartPoint==0) {
                            generatorStartPoint=time(0);
                            printf("Setting start position of p-random number generator to actual CPU time...\n");
                        } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                        InitMT((unsigned int)generatorStartPoint);
                        randomStartStep[0]=generatorStartPoint;
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
                printf("Start of equilibration at reduced pressure: %.7E (startDen: %.7E, startPacFrac: %.7E)...",pressureReduced,rho,pacFrac);
                if (autoEqLength==0) printf(" (%ld cycles)\n",cyclesOfEquilibration); else printf(" (auto, min. length: %ld cycles)\n",cyclesOfEquilibration*intervalSampling);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.7E\n",N,gaps,growing,pressureReduced);
            fflush(stdout);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                fullCycle=activeN+volumeMoveChance,fullCycle2=fullCycle*2,fullCycleDiagonalMove=fullCycle2-volumeMoveChance,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberPhi=0, displacedNumberPhi=0,
                attemptedNumberV=0, displacedNumberV=0,
                attemptedNumberVND=0, displacedNumberVND=0,
                cyclesOfMeasurementBuffer=arg[1]==0?cyclesOfMeasurement:0,
                cyclesOfEquilibrationBuffer=cyclesOfEquilibration;
            double deltaRTable[100], deltaRMean=deltaR, deltaPhiTable[100], deltaPhiMean=deltaPhi, deltaVTable[100], deltaVMean=deltaV, deltaVNDTable[100], deltaVNDMean=deltaVND;
            for (int i=0;i<100;i++) {deltaRTable[i]=deltaRMean; deltaPhiTable[i]=deltaPhiMean; deltaVTable[i]=deltaVMean; deltaVNDTable[i]=deltaVNDMean;}
            int simulationStage=cyclesOfEquilibrationBuffer>0?0:cyclesOfMeasurementBuffer>0?1:2;  //0-equilibration, 1-measurement, 2-end
            double autoEqVolumeTable[cyclesOfEquilibrationBuffer],autoEqBalance=0;
            autoEqVolumeTable[0]=boxMatrixInnerParameters[1]; for (int i=1;i<cyclesOfEquilibrationBuffer;i++) autoEqVolumeTable[i]=-1;
            int cycleCounter=0, indexScanned=(matrixOfParticlesSize[0]*round(matrixOfParticlesSize[1]*NLinearMod/2.0)-round(matrixOfParticlesSize[0]/2.0))*NLinearMod,autoEqCounter=0,autoEqCheck=0;  //TODO: indexScanned w 3D inaczej powinien być ustawiany
            bool volumeMove=false;

            char allResultsFileName[200],bufferConfigurations[200],bufferSavedConfigurations[200],bufferOrientations[200],bufferOrientatCorrelFun[200],allOrientationsFileName[200],bufferOrientationsResults[200],allOrientationsResultsFileName[200],bufferPressure[100];
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferOrientations,orientationsFileName); strcpy(allOrientationsFileName,orientationsFileName);
            if (OCFMode) strcpy(bufferOrientatCorrelFun,orientatCorrelFunFileName);
            strcpy(bufferOrientationsResults,orientationsResultsFileName); strcpy(allOrientationsResultsFileName,orientationsResultsFileName);
            sprintf(bufferPressure,"%.4E",pressureReduced);
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
            if (onlyMath[0]) simulationStage=2;

            timeStart=time(0);
            while (simulationStage<2) {
                int randIndex;
                if (volumeMove) {
                    randIndex = (int)(MTRandom0to1(randomStartStep)*activeN2);
                    volumeMove=false;
                } else randIndex = (int)(MTRandom0to1(randomStartStep)*fullCycle2);
                if (randIndex<activeN) {
                    attemptedNumberR++;
                    if (attemptToDisplaceAParticle(atoms,particles,randIndex,boxMatrix,boxMatrixInnerParameters,&totalInteractionEnergy,true))
                        displacedNumberR++;
                } else if (randIndex<activeN2) {
                    attemptedNumberPhi++;
                    if (attemptToDisplaceAParticle(atoms,particles,randIndex-activeN,boxMatrix,boxMatrixInnerParameters,&totalInteractionEnergy,false))
                        displacedNumberPhi++;
                } else {
                    volumeMove=true;
                    if (randIndex<fullCycleDiagonalMove) {
                        attemptedNumberV++;
                        if (attemptToChangeVolumeSeparate(atoms,particles,pressure,boxMatrix,boxMatrixInnerParameters,&totalInteractionEnergy,deltaV,true))
                            displacedNumberV++;
                    } else {
                        attemptedNumberVND++;
                        if (attemptToChangeVolumeSeparate(atoms,particles,pressure,boxMatrix,boxMatrixInnerParameters,&totalInteractionEnergy,deltaVND,false))
                            displacedNumberVND++;
                    }
                }

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++;

                    if (cycle%intervalSampling==0) {
                        if (simulationStage==0) {
                            if (autoEqLength==0) {
                                if (cycle>cyclesOfEquilibrationBuffer) simulationStage=1;
                            } else {  //metoda 'przyrostow' w sposob 'listy' przy uzyciu tablicy-ale bez przepisywania wszystkich elementow (o 1 w dol) jezeli tablica sie zapelnia
                                autoEqBalance+=boxMatrixInnerParameters[1]-autoEqVolumeTable[autoEqCounter++];
                                if (autoEqCounter>=cyclesOfEquilibrationBuffer) autoEqCounter=0;
                                int nextIndex=autoEqCounter+1; if (nextIndex>=cyclesOfEquilibrationBuffer) {nextIndex=0; if (autoEqCheck==0) autoEqCheck=autoEqBalance>0?1:-1;}
                                if (autoEqVolumeTable[autoEqCounter]>0) autoEqBalance-=autoEqVolumeTable[nextIndex]-autoEqVolumeTable[autoEqCounter];
                                autoEqVolumeTable[autoEqCounter]=boxMatrixInnerParameters[1];
                                if (autoEqCheck!=0 && autoEqBalance*autoEqCheck<0) {
                                    simulationStage=1; cyclesOfEquilibrationBuffer=cycle;
                                }
                            }
                            if (simulationStage==1) {
                                printf("Equilibration finished after: %ld cycles (%ldsec).\n",cyclesOfEquilibrationBuffer,time(0)-timeStart);
                                fflush(stdout);
                            }
                        }
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR,
                               acceptanceRatioPhi = displacedNumberPhi/(double)attemptedNumberPhi,
                               acceptanceRatioV = displacedNumberV/(double)attemptedNumberV,
                               acceptanceRatioVND = displacedNumberVND/(double)attemptedNumberVND;
                        if (cycle%neighUpdatingFrequency==0) {
                            updateAtomsNeighbourList(atoms,particles,boxMatrix,boxMatrixInnerParameters[1],true);
                            totalInteractionEnergy=getEnergyAll(atoms,particles,boxMatrix);
                        }

                        /////wypisywanie danych czesciej niz normalnie i PRZED zrownowagowaniem
                        /*if (cycle%50==0) {
                            rho=N2/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                            printf("Cycle: %ld\n",(cycle+args[4]));
                            printf("   AccRatR: %.4E, dR: %.4E, AccRatPhi: %.4E, dPhi: %.4E\n   AccRatV: %.4E, dV: %.4E, AccRatVND: %.4E, dVND: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioPhi,deltaPhi,acceptanceRatioV,deltaV,acceptanceRatioVND,deltaVND);
                            printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                            printf("   box00: %.8E, box11: %.8E, box22: %.8E\n   box01(10): %.8E, box02(20): %.8E, box12(21): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                            fflush(stdout);
                        }*/
                        /////

                        if (simulationStage==1) {
                            if (timeEquilibration==0) timeEquilibration=time(0);

                            rho=N2/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0)
                                fprintf(fileAllResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);

                            if (cycle%intervalOrientations==0) {
                                /////skan po konkretnej cząstce (np. dla 224: 14*(16/2)-(14/2)=105, etc.) - w srodku by nie skakala na granicy pudla periodycznego; bezposrednie uzycie w Mathematice (format tablicy)
                                if (cycle-cyclesOfEquilibrationBuffer>=cyclesOfMeasurementBuffer)
                                    fprintf(fileOrientations,"{%.12E,%.12E,%.12E,%.12E,%.12E,%.12E}}",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].r[2],particles[indexScanned].phi[0],particles[indexScanned].phi[1],particles[indexScanned].phi[2]);
                                else fprintf(fileOrientations,"{%.12E,%.12E,%.12E,%.12E,%.12E,%.12E},",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].r[2],particles[indexScanned].phi[0],particles[indexScanned].phi[1],particles[indexScanned].phi[2]);
                                /////skan po wszystkich cząstkach
                                for (int i=0;i<activeN-1;i++) fprintf(fileAllOrientations,"{%.12E,%.12E,%.12E},",particles[i].phi[0],particles[i].phi[1],particles[i].phi[2]);
                                fprintf(fileAllOrientations,"{%.12E,%.12E,%.12E}\n",particles[activeN-1].phi[0],particles[activeN-1].phi[1],particles[activeN-1].phi[2]);
                                /////OCF
                                /*if (OCFMode) {
                                    char bufferText[4096]="", bufferAngle[20];
                                    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                        if (i<particles[i].neighbours[j]) {
                                            strcpy(bufferText,"");
                                            getParticlesDistanceSquared(particles[i],particles[particles[i].neighbours[j]],boxMatrix);
                                            double gamma=atan(dr[1]/dr[0]),
                                                   aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                                                   bAngle=particles[i].phi-gamma;
                                            aAngle=normalizeAngle(aAngle+C); bAngle=normalizeAngle(bAngle+C);

                                            strncat(bufferText,"{",2); sprintf(bufferAngle,"%.12E",aAngle); strncat(bufferText,bufferAngle,20);
                                            strncat(bufferText,",",2); sprintf(bufferAngle,"%.12E",bAngle); strncat(bufferText,bufferAngle,20);
                                            if (i>=activeN-2 && cycle-cyclesOfEquilibrationBuffer>=cyclesOfMeasurementBuffer) strncat(bufferText,"}}",3); // -2 because it's last particle which has neighbour UNtested (the last one)
                                            else strncat(bufferText,"},",3);
                                            fprintf(fileOrientatCorrelFun,"%s",bufferText);
                                        }
                                    }
                                }*/
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t{",(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++) fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E,%.17E,%.17E,%.17E],",particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].phi[0],particles[i].phi[1],particles[i].phi[2]);
                                fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E,%.17E,%.17E,%.17E]}",particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].r[2],particles[activeNMinus1].phi[0],particles[activeNMinus1].phi[1],particles[activeNMinus1].phi[2]);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                printf("Cycle: %ld, ",(cycle+(long)args[4]));
                                printf("simulation time: full-%ldsec, measurement-%ldsec\n",time(0)-timeStart,time(0)-timeEquilibration);
                                printf("   AccRatR: %.4E, dR: %.4E, AccRatPhi: %.4E, dPhi: %.4E\n   AccRatV: %.4E, dV: %.4E, AccRatVND: %.4E, dVND: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioPhi,deltaPhi,acceptanceRatioV,deltaV,acceptanceRatioVND,deltaVND);
                                printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8E, box11: %.8E, box22: %.8E\n   box01(10): %.8E, box02(20): %.8E, box12(21): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                                printf("   totalInterEnergy: %.8E\n",totalInteractionEnergy);
                                fflush(stdout);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2],deltaR,deltaPhi,deltaV,deltaVND);
                                for (int i=0;i<activeN;i++) fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.17E,%.17E,%.17E,%.17E,%.17E,%.17E],",particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].phi[0],particles[i].phi[1],particles[i].phi[2],particles[i].halfSigma,particles[i].atoms[0]->diameter,particles[i].atoms[1]->diameter);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Parallelepiped[{0,0,0},{{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E}}]},{Opacity[0.2],Green,Sphere[{%.12E,%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][0],boxMatrix[1][1],boxMatrix[1][2],boxMatrix[2][0],boxMatrix[2][1],boxMatrix[2][2],particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].r[2],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[2][2]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E, boxMatrix[0][2]=boxMatrix[2][0]=%.12E, boxMatrix[1][2]=boxMatrix[2][1]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                                fclose(fileConfigurations);
                            }
                        } else {//dostosowywanie delt w trakcie pomiaru
                            if (acceptanceRatioR>desiredAcceptanceRatioR) {
                                deltaR*=1.05;
                                if (deltaR>maxDeltaR) deltaR=maxDeltaR;  //jednostka odległości to 1; ograniczenie od góry, żeby cząstka nie skakała zbyt mocno po pudle, nie teleportowała się
                            } else {
                                deltaR*=0.95;
                                while (boxMatrix[0][0]+deltaR*0.1==boxMatrix[0][0]) deltaR*=1.05;  //dolna granica deltaR - gdy jest tak małe (w stosunku do boxMatrix[0][0]), że go nie może zmienić ze względu na ograniczoną precyzję (czyli, nie może zmienić pozycji cząstki znajdujące1j się przy granicy pudła), należy ją zwiększać tak długo, aż będzie mogła wprowadzać takie zmiany. Bez uwzględnienia dolnej granicy delt, po jej przekroczeniu, dana delta przestaje być wrażliwa na dostosowywanie (czy się ją delikatnie zwiększy czy zmniejszy, nie zmieni to poziomu akceptacji). Acceptance ruchów będzie wówczas zazwyczaj dość mały (nie 50%, znacznie mniejszy, mniejszy niż 20% - jeżeli byłby duży, to nigdy by do takiej sytuacji nie doszło - bo przecież delta byłaby zwiększana a nie zmniejszana) i będzie dochodziło do dalszego zmniejszania danej delty, aż do -INF. Mnożenie *0.1: 0.5(zakres losowy od -0.5 do 0.5)*0.5(średnia wartość zmiennej losowej od 0 do 1)*0.4(tak dla bezpieczeństwa, żeby nie schodzić do granicznie małych wartości)
                            }
                            if (acceptanceRatioPhi>desiredAcceptanceRatioR) {
                                deltaPhi*=1.05;                         //ograniczenie od góry, aby dimer nie mógł dokonać 'skokowego obrotu' w gęstym upakowaniu, co sprawiałoby, że atomy zamieniałyby sąsiadów
                                if (deltaPhi>maxDeltaR) deltaPhi=maxDeltaR;
                            } else {
                                deltaPhi*=0.95;
                                while (boxMatrix[0][0]+deltaPhi*0.1==boxMatrix[0][0]) deltaPhi*=1.05;   //ograniczenie od dołu jak dla deltaR - trochę nieuzasadnione, ale deltaR i deltaPhi przyjmują bardzo podobne wartości, można zatem przyjąć, że taka sama blokada będzie odpowiednia
                            }
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; //ograniczenie deltaV od góry chyba nie jest potrzebne, w płynie kształt i tak jest blokowany, same ruchy objętościowe nie wprowadzą niczego złego, nawet jak będą bardzo duże (dalej są określane przez entalpię swobodną)
                            else {
                                deltaV*=0.95;
                                while (exp(log(boxMatrix[0][0])+deltaV*0.1)==boxMatrix[0][0]) deltaV*=1.05;     //ograniczenie od dołu jak dla deltaR (ruchy logarytmiczne elementów diagonalnych)
                            }
                            if (acceptanceRatioVND>desiredAcceptanceRatioV) deltaVND*=1.05;
                            else {
                                deltaVND*=0.95;
                                while (boxMatrix[0][1]+deltaVND*0.1==boxMatrix[0][1]) deltaVND*=1.05;   //ograniczenie od dołu jak dla deltaR, ale bazujące na przykładowym elemencie niediagonalnym [tu: 01] (ruchy liniowe elementów niediagonalnych)
                            }

                            int sampleNumberMod100=(cycle/intervalSampling)%100;
                            updateTableAndGetActualMean(deltaRTable,deltaRMean,sampleNumberMod100,deltaR); deltaR=deltaRMean;
                            updateTableAndGetActualMean(deltaPhiTable,deltaPhiMean,sampleNumberMod100,deltaPhi); deltaPhi=deltaPhiMean;
                            updateTableAndGetActualMean(deltaVTable,deltaVMean,sampleNumberMod100,deltaV); deltaV=deltaVMean;
                            updateTableAndGetActualMean(deltaVNDTable,deltaVNDMean,sampleNumberMod100,deltaVND); deltaVND=deltaVNDMean;
                        }
                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberPhi=0; displacedNumberPhi=0;
                        attemptedNumberV=0; displacedNumberV=0;
                        attemptedNumberVND=0; displacedNumberVND=0;
                    }
                    if (simulationStage==1 && cycle-cyclesOfEquilibrationBuffer>=cyclesOfMeasurementBuffer) simulationStage=2;
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
            char linia[4096]; double dataLicznik=0; int faultyLines=0, onlyMathLinesBuffer=onlyMath[1];
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) {
                if (fgets(linia,300,fileAllResults)!=NULL && !isLineCorrect(linia)) onlyMathLinesBuffer++;
            }
            while(fgets(linia,300,fileAllResults)!=NULL) {
                if (isLineCorrect(linia)) dataLicznik++; else faultyLines++;
            }
            fclose(fileAllResults);
            printf("done (Found %ld data lines [%d faulty lines occurred].",(long)dataLicznik,faultyLines);
            if ((long)dataLicznik%10>0) printf(" Last %ld lines won't be considered, due to calculations of averages in 10 sets.)\n",(long)dataLicznik%10); else printf(")\n");
            dataLicznik-=(long)dataLicznik%10;

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... "); fflush(stdout);
            double avVolumeSet[10], avBoxMatrixSet[10][6], avRhoSet[10], avPacFracSet[10]; //wyniki dzielone są na 10 zespołów (obliczane są nieskorelowane wzajemnie średnie "lokalne", do obliczenia błędu średniej "globalnej")
            for (int i=0;i<10;i++) {
                avVolumeSet[i]=0; for (int j=0;j<6;j++) avBoxMatrixSet[i][j]=0;
                avRhoSet[i]=0; avPacFracSet[i]=0;
            }
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            double lineCounter=0;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                int dataIndex=0; double dataD[6]; while (dataIndex<6) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10;
                    else dataD[dataIndex++]=strtod(data,NULL);
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    for (int i=0;i<6;i++) avBoxMatrixSet[setIndex][i]+=dataD[i];
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            double avVolume=0, avBoxMatrix[6]={0,0,0,0,0,0}, avRho=0, avPacFrac=0;
            for (int i=0;i<10;i++) {
                for (int j=0;j<6;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
                avVolumeSet[i]=fabs(avBoxMatrixSet[i][0]*avBoxMatrixSet[i][1]*avBoxMatrixSet[i][2]+avBoxMatrixSet[i][3]*avBoxMatrixSet[i][5]*avBoxMatrixSet[i][4]+
                                    avBoxMatrixSet[i][4]*avBoxMatrixSet[i][3]*avBoxMatrixSet[i][5]-avBoxMatrixSet[i][4]*avBoxMatrixSet[i][1]*avBoxMatrixSet[i][4]-
                                    avBoxMatrixSet[i][0]*avBoxMatrixSet[i][5]*avBoxMatrixSet[i][5]-avBoxMatrixSet[i][3]*avBoxMatrixSet[i][3]*avBoxMatrixSet[i][2]);
                avVolume+=avVolumeSet[i];
                avRhoSet[i]=N2/avVolumeSet[i]; avRho+=avRhoSet[i];
                avPacFracSet[i]=1.0/VcpPerParticle/avRhoSet[i]; avPacFrac+=avPacFracSet[i];
            }
            avVolume*=0.1; avRho*=0.1; avPacFrac*=0.1; for (int i=0;i<6;i++) avBoxMatrix[i]*=0.1;
            //obliczenie bledow mierzonych wielkosci
            double dAvVolume=0, dAvBoxMatrix[6]={0,0,0,0,0,0}, dAvRho=0, dAvPacFrac=0;
            for (int i=0;i<10;i++) {double epsilon=avVolume-avVolumeSet[i]; dAvVolume+=epsilon*epsilon;} dAvVolume=getAvErrorFromSumEps(dAvVolume,90.0); //10*9 (n(n-1))
            for (int j=0;j<6;j++) {for (int i=0;i<10;i++) {double epsilon=avBoxMatrix[j]-avBoxMatrixSet[i][j]; dAvBoxMatrix[j]+=epsilon*epsilon;} dAvBoxMatrix[j]=getAvErrorFromSumEps(dAvBoxMatrix[j],90.0);}
            for (int i=0;i<10;i++) {double epsilon=avRho-avRhoSet[i]; dAvRho+=epsilon*epsilon;} dAvRho=getAvErrorFromSumEps(dAvRho,90.0);
            for (int i=0;i<10;i++) {double epsilon=avPacFrac-avPacFracSet[i]; dAvPacFrac+=epsilon*epsilon;} dAvPacFrac=getAvErrorFromSumEps(dAvPacFrac,90.0);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... "); fflush(stdout);
            double e0000Set[10], e0001Set[10], e0002Set[10], e0011Set[10], e0012Set[10], e0022Set[10],
                   e0101Set[10], e0102Set[10], e0111Set[10], e0112Set[10], e0122Set[10],
                   e0202Set[10], e0211Set[10], e0212Set[10], e0222Set[10],
                   e1111Set[10], e1112Set[10], e1122Set[10],
                   e1212Set[10], e1222Set[10],
                   e2222Set[10],
                   H00,H11,H22,H01,H02,H12,
                   H00_2,H11_2,H22_2,H01_2,H01_3,H01_4,H02_2,H02_3,H02_4,H12_2,H12_3,H12_4,
                   H12_2mH11tH22,H11tH12,H00tH12,H00tH11,H01tH12,H01tH02,H02tH11,H02tH12,H11tH22,H00tH02,
                   denominatorComponent,denominatorInv;
            for (int i=0;i<10;i++) {
                e0000Set[i]=0; e0001Set[i]=0; e0002Set[i]=0; e0011Set[i]=0; e0012Set[i]=0; e0022Set[i]=0;
                e0101Set[i]=0; e0102Set[i]=0; e0111Set[i]=0; e0112Set[i]=0; e0122Set[i]=0;
                e0202Set[i]=0; e0211Set[i]=0; e0212Set[i]=0; e0222Set[i]=0;
                e1111Set[i]=0; e1112Set[i]=0; e1122Set[i]=0;
                e1212Set[i]=0; e1222Set[i]=0;
                e2222Set[i]=0;
            }
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            lineCounter=0; int setIndex, oldSetIndex=-1;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                setIndex=(int)(lineCounter/dataLicznik*10.0);
                if (setIndex!=oldSetIndex) {
                    H00=avBoxMatrixSet[setIndex][0]; H11=avBoxMatrixSet[setIndex][1]; H22=avBoxMatrixSet[setIndex][2];
                    H01=avBoxMatrixSet[setIndex][3]; H02=avBoxMatrixSet[setIndex][4]; H12=avBoxMatrixSet[setIndex][5];
                    H00_2=H00*H00; H11_2=H11*H11; H22_2=H22*H22; H01_2=H01*H01; H01_3=H01*H01_2; H01_4=H01_2*H01_2; H02_2=H02*H02; H02_3=H02*H02_2; H02_4=H02_2*H02_2; H12_2=H12*H12; H12_3=H12*H12_2; H12_4=H12_2*H12_2;
                    H12_2mH11tH22=H12_2-H11*H22; H11tH12=H11*H12; H00tH12=H00*H12; H00tH11=H00*H11; H01tH12=H01*H12; H01tH02=H01*H02; H02tH11=H02*H11; H02tH12=H02*H12; H11tH22=H11*H22; H00tH02=H00*H02;
                    denominatorComponent=H02_2*H11-2*H01tH02*H12+H01_2*H22+H00*H12_2mH11tH22; denominatorInv=1/(2*denominatorComponent*denominatorComponent);
                    oldSetIndex=setIndex;
                }
                sscanf(linia,"%c",linia);
                int actIndex=0;
                double h[6],h00,h11,h22,h01,h02,h12;
                int dataIndex=0; while (dataIndex<6) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10) {
                        if ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.')) dataIndex=10;
                        else h[dataIndex++]=strtod(data,NULL);
                    }
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    h00=h[0]; h11=h[1]; h22=h[2]; h01=h[3]; h02=h[4]; h12=h[5];
                    double h00_2=h00*h00,h11_2=h11*h11,h22_2=h22*h22,h01_2=h01*h01,h02_2=h02*h02,h12_2=h12*h12,
                           h12th22=h12*h22,h01th02=h01*h02,h00th01=h00*h01,h11th12=h11*h12,h01th12=h01*h12,h02th12=h02*h12,h00th02=h00*h02,h01th11=h01*h11,
                           h00ph11=h00+h11,h00ph22=h00+h22,h11ph22=h11+h22,

                           e00=(-H02_4*H11_2+4*H01*H02_3*H11tH12+H12_2*(-2*h01*H01*h12*H12+(h00_2-H00_2+h01_2)*H12_2+H01_2*(h12_2+h22_2))-2*H12*((h00_2-H00_2+h01_2)*H11tH12-h01*H01*(H11*h12+(h00ph11)*H12)+H01_2*(H00tH12+h12*(h11ph22)))*H22+(-H01_4-2*h01*H01*(h00ph11)*H11+(h00-H00)*(h00+H00)*H11_2+h01_2*(H01_2+H11_2)+H01_2*(h11_2+2*H00tH11+h12_2))*H22_2+H02_2*((h01_2-4*H01_2+h11_2+h12_2)*H12_2+H11_2*(h12_2+h22_2+2*H00*H22)-2*H11*(H12*(H00tH12+h12*(h11ph22))+H01_2*H22))+2*H02*(H01tH12*(-H11*(h12_2+h22_2)+H12*(2*H00tH12+h12*(h11ph22)))+2*H01_3*H12*H22-H01*(-h11*H11*h12+h11_2*H12+(h01_2+2*H00tH11+h12_2)*H12-H11*h12th22)*H22+h01*(H11*h12-(h00ph11)*H12)*H12_2mH11tH22)+h02_2*(H02_2*H11_2-2*H01tH02*H11tH12+H01_2*H12_2+H12_2mH11tH22*H12_2mH11tH22)+2*h02*(-h01*(H02tH11-H01tH12)*(H02tH12-H01*H22)-H12_2mH11tH22*(-h00*H02tH11+h00*H01tH12+H02*h12*H12-H02tH11*h22+H01tH12*h22-H01*h12*H22)))*denominatorInv,
                           e11=(h01_2*H02_4+H02_4*h11_2-H02_4*H11_2+H02_4*h12_2+2*H00*h01th02*H02_2*H12-2*h00th01*H02_3*H12-2*h01*H02_3*h11*H12-2*h02*H02_3*h12*H12+2*H00*H02_2*h11th12*H12+H00_2*h02_2*H12_2-2*h00*H00*h02*H02*H12_2+h00_2*H02_2*H12_2+h01_2*H02_2*H12_2+h02_2*H02_2*H12_2-2*H00*H02_2*H11*H12_2-2*H00*h01*H02*h12*H12_2+H00_2*h12_2*H12_2-H00_2*H12_4+2*H00*H02_2*h12*H12*h22-2*H00*h02*H02*H12_2*h22+H00_2*H12_2*h22_2+4*H01_3*H02tH12*H22-2*H00*(h01_2*H02_2+H02_2*(h11_2-H11_2+h12_2)+h01*(H00*h02-H02*(h00ph11))*H12-h02*H02*h12*H12+H00tH12*(-H11tH12+h12*(h11ph22)))*H22-H01_4*H22_2+H00_2*(h01_2+h11_2-H11_2+h12_2)*H22_2+H01_2*(H02_2*(h02_2+h12_2-4*H12_2+h22_2)-2*(H00*H12_2+H02*(H02tH11+h01th12+h02*(h00ph22)))*H22+(h00_2+h01_2+h02_2+2*H00tH11)*H22_2)+2*H01*(-H02*(h01*H02*(h02*H02-h12*H12)-h02*H02tH12*(h00ph22)+H00tH12*(h02_2+h12_2-2*H12_2+h22_2)+H02_2*(-2*H11tH12+h12*(h11ph22)))+(H02*(h01*H02*(h00ph11)+h02*H02*h12-(h00_2+h01_2+h02_2)*H12)+H00*(h01*(h02*H02+h12*H12)+h02*H12*(h00ph22)+H02*(-2*H11tH12+h12*(h11ph22))))*H22-H00*(h01*(h00ph11)+h02th12)*H22_2))*denominatorInv,
                           e22=(H02_2*(h00_2+h01_2+h02_2-H02_2)*H11_2+2*H00*H02tH11*(-h00th02*H11-h01*H11*h12+h00th01*H12+h01th11*H12+h02th12*H12-H02*H12_2-h02*H11*h22+H02tH11*H22)-2*H01_3*(h01*(h02*H02+h12*H12)+h02*H12*(h00ph22)+H02*h12*(h11ph22)-2*H02tH12*H22)+H01_4*(h02_2+h12_2+h22_2-H22_2)+H01_2*(-2*H00*h02_2*H11-2*H00tH11*h12_2+2*H00*h11th12*H12+h00_2*H12_2+h02_2*H12_2+2*h01*(H02tH11*h12+H00*h02*H12+H02*(h00ph11)*H12)+h01_2*(H02_2+H12_2)+2*H00*h12*H12*h22-2*H00tH11*h22_2+2*h02*H02*(h12*H12+H11*(h00ph22))-2*H00*H12_2*H22+2*H00tH11*H22_2+H02_2*(h11_2+h12_2-4*H12_2-2*H11tH22))+H00_2*(h02_2*H11_2-2*h01th02*H11tH12+H12_2*(h01_2+h11_2+h12_2-H12_2)+2*H11tH12*(-h12*(h11ph22)+H12*H22)+H11_2*(h12_2+h22_2-H22_2))-2*H01*(H02tH11*(h01*H02*(h00ph11)+h02*H02*h12+(h00_2+h01_2+h02_2-2*H02_2)*H12)+H00*(h01_2*H02tH12+h01*(-h02*H02tH11+H12*(-H11*h12+(h00ph11)*H12))+h02*H12*(h12*H12-H11*(h00ph22))+H02*(-h11*H11*h12+h11_2*H12+h12_2*H12-2*H12_3-H11*h12th22+2*H11tH12*H22))))*denominatorInv,
                           e01=(h01th02*H02_3*H11+H02_3*h11*H11*h12-h01_2*H02_3*H12-H02_3*h11_2*H12+H00*h02_2*H02tH11*H12-h00th02*H02_2*H11tH12-h01*H02_2*H11*h12*H12-H02_3*h12_2*H12+H00*H02tH11*h12_2*H12-H00*h01th02*H02*H12_2+2*h00th01*H02_2*H12_2+2*h01*H02_2*h11*H12_2+2*h02*H02_2*h12*H12_2-H00tH02*h11th12*H12_2+h00*H00*h02*H12_3-h00_2*H02*H12_3-h01_2*H02*H12_3-h02_2*H02*H12_3+H00*h01th12*H12_3+H02_3*H11*h12th22-h02*H02_2*H11tH12*h22-H00tH02*h12*H12_2*h22+H00*h02*H12_3*h22+H00*H02tH11*H12*h22_2+(H02tH11*(-H02*(h01*(h00ph11)+h02th12)+(h00_2+h01_2+h02_2)*H12)+H00*(h01_2*H02tH12-h01*(h02*H02tH11+H12*(H11*h12+(h00ph11)*H12))+H02*(-h11*H11*h12+h11_2*H12+h12_2*H12-H11*h12th22)-h02*H12*(h12*H12+H11*(h00ph22))))*H22+H00tH11*(h01*(h00ph11)+h02th12)*H22_2+H01_2*(H02tH12*(h02_2+h12_2+h22_2)-(h01*(h02*H02+h12*H12)+h02*H12*(h00ph22)+H02*h12*(h11ph22))*H22+(h01*(h00ph11)+h02th12)*H22_2)+H01*(-(H02_2*H11+H00*H12_2)*(h02_2+h12_2+h22_2)+(H02_2*(h11_2+h12_2)+2*h01*(H02tH11*h12+H00*h02*H12-H02*(h00ph11)*H12)+h01_2*(H02_2+H12_2)+2*h02*H02*(-h12*H12+H11*(h00ph22))+H12*((h00_2+h02_2)*H12+2*H00*h12*(h11ph22)))*H22-((h00_2+h01_2+h02_2)*H11+H00*(h01_2+h11_2+h12_2))*H22_2))*denominatorInv,
                           e02=(-H00*(h02_2*H02*H11_2+h01*H12_2*(H11*h12-(h00ph11)*H12)+H02*((h01_2+h11_2+h12_2)*H12_2-2*H11*h12*H12*(h11ph22)+H11_2*(h12_2+h22_2))+h02*H12*(-2*h01*H02tH11+H12*(-h12*H12+H11*(h00ph22))))+H00tH11*(h00th02*H11+h01*H11*h12-h00th01*H12-h01th11*H12-h02th12*H12+h02*H11*h22)*H22+H02tH11*(h00th02*H02tH11+h01*H02tH11*h12-h00th01*H02tH12-h01*H02*h11*H12-h02*H02*h12*H12+h00_2*H12_2+h01_2*H12_2+h02_2*H12_2+h02*H02tH11*h22-(h00_2+h01_2+h02_2)*H11tH22)+H01_2*(h02_2*H02tH11+2*h01th12*H12_2+2*h02*H12_2*(h00ph22)+H02tH11*(h12_2+h22_2)-H02*(h01_2+h11_2+h12_2)*H22-h01*(H11*h12+(h00ph11)*H12)*H22-h02*(h12*H12+H11*(h00ph22))*H22)+H01_3*(-H12*(h02_2+h12_2+h22_2)+(h01th02+h12*(h11ph22))*H22)+H01*(H02_2*(-h11*H11*h12+h11_2*H12+h12*(h12*H12-H11*h22))-H12*((h00_2+h02_2)*H12_2-H00*(h02_2*H11-h12*H12*(h11ph22)+H11*(h12_2+h22_2)))+((h00_2+h02_2)*H11tH12+H00*(-h11*H11*h12+h11_2*H12+h12*(h12*H12-H11*h22)))*H22+h01_2*H12*(H02_2-H12_2+(H00+H11)*H22)-2*h02*H02tH11*(H12*(h00ph22)-h12*H22)-h01*(2*H02tH11*(h12*H12-(h00ph11)*H22)+h02*(H02_2*H11+H00*(H12_2+H11tH22)))))*denominatorInv,
                           e12=(H02_2*H11*(h01*H02*(h00ph11)+h02*H02*h12-(h00_2+h01_2+h02_2)*H12)+H01_3*(-H02*(h02_2+h12_2+h22_2)+(h01th12+h02*(h00ph22))*H22)+H01_2*(2*h01th02*H02_2+2*H02_2*h11th12+H00*h02_2*H12+H00*h12_2*H12+2*H02_2*h12th22+H00tH12*h22_2-(h01*H02*(h00ph11)+h02*H02*h12+(h00_2+h01_2+h02_2)*H12+H00*(h01th02+h12*(h11ph22)))*H22)+H00tH02*(h01_2*H02tH12+h02*H12*(-h12*H12+2*H11*(h00ph22))+H02*(-h11*H11*h12+h11_2*H12+h12*(h12*H12-H11*h22))-h02*H11*h12*H22-h01*(h02*H02tH11-2*H11*h12*H12+(h00ph11)*H12_2+(h00ph11)*H11tH22))+H00_2*(-h02_2*H11tH12-H11tH12*(h12_2+h22_2)+H11*h12*(h11ph22)*H22+h01th02*(H12_2+H11tH22)+H12*(h12*H12*(h11ph22)-(h01_2+h11_2+h12_2)*H22))+H01*(-H02_3*(h11_2+h12_2)-h02*H02_2*H11*(h00ph22)+h01_2*H02*(-H02_2+H12_2+(H00+H11)*H22)-H00*h02*(H12_2*(h00ph22)+(-2*h12*H12+H11*(h00ph22))*H22)+H02*((h00_2+h02_2)*(H12_2+H11tH22)+H00*(h02_2*H11-2*h12*H12*(h11ph22)+H11*(h12_2+h22_2)+(h11_2+h12_2)*H22))-h01*(H02_2*H11*h12+2*H00*h02*H02tH12+H00*(-2*(h00ph11)*H12*H22+h12*(H12_2+H11tH22)))))*denominatorInv;

                    e0000Set[setIndex]+=e00*e00; e0001Set[setIndex]+=e00*e01; e0002Set[setIndex]+=e00*e02; e0011Set[setIndex]+=e00*e11; e0012Set[setIndex]+=e00*e12; e0022Set[setIndex]+=e00*e22;
                    e0101Set[setIndex]+=e01*e01; e0102Set[setIndex]+=e01*e02; e0111Set[setIndex]+=e01*e11; e0112Set[setIndex]+=e01*e12; e0122Set[setIndex]+=e01*e22;
                    e0202Set[setIndex]+=e02*e02; e0211Set[setIndex]+=e02*e11; e0212Set[setIndex]+=e02*e12; e0222Set[setIndex]+=e02*e22;
                    e1111Set[setIndex]+=e11*e11; e1112Set[setIndex]+=e11*e12; e1122Set[setIndex]+=e11*e22;
                    e1212Set[setIndex]+=e12*e12; e1222Set[setIndex]+=e12*e22;
                    e2222Set[setIndex]+=e22*e22;
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            double e0000=0, e0001=0, e0002=0, e0011=0, e0012=0, e0022=0,
                   e0101=0, e0102=0, e0111=0, e0112=0, e0122=0,
                   e0202=0, e0211=0, e0212=0, e0222=0,
                   e1111=0, e1112=0, e1122=0,
                   e1212=0, e1222=0,
                   e2222=0;
            for (int i=0;i<10;i++) {
                e0000Set[i]/=dataLicznik*0.1; e0000+=e0000Set[i]; e0001Set[i]/=dataLicznik*0.1; e0001+=e0001Set[i]; e0002Set[i]/=dataLicznik*0.1; e0002+=e0002Set[i];
                e0011Set[i]/=dataLicznik*0.1; e0011+=e0011Set[i]; e0012Set[i]/=dataLicznik*0.1; e0012+=e0012Set[i]; e0022Set[i]/=dataLicznik*0.1; e0022+=e0022Set[i];
                e0101Set[i]/=dataLicznik*0.1; e0101+=e0101Set[i]; e0102Set[i]/=dataLicznik*0.1; e0102+=e0102Set[i]; e0111Set[i]/=dataLicznik*0.1; e0111+=e0111Set[i];
                e0112Set[i]/=dataLicznik*0.1; e0112+=e0112Set[i]; e0122Set[i]/=dataLicznik*0.1; e0122+=e0122Set[i];
                e0202Set[i]/=dataLicznik*0.1; e0202+=e0202Set[i]; e0211Set[i]/=dataLicznik*0.1; e0211+=e0211Set[i]; e0212Set[i]/=dataLicznik*0.1; e0212+=e0212Set[i];
                e0222Set[i]/=dataLicznik*0.1; e0222+=e0222Set[i];
                e1111Set[i]/=dataLicznik*0.1; e1111+=e1111Set[i]; e1112Set[i]/=dataLicznik*0.1; e1112+=e1112Set[i]; e1122Set[i]/=dataLicznik*0.1; e1122+=e1122Set[i];
                e1212Set[i]/=dataLicznik*0.1; e1212+=e1212Set[i]; e1222Set[i]/=dataLicznik*0.1; e1222+=e1222Set[i];
                e2222Set[i]/=dataLicznik*0.1; e2222+=e2222Set[i];
            }
            e0000*=0.1; e0001*=0.1; e0002*=0.1; e0011*=0.1; e0012*=0.1; e0022*=0.1;
            e0101*=0.1; e0102*=0.1; e0111*=0.1; e0112*=0.1; e0122*=0.1;
            e0202*=0.1; e0211*=0.1; e0212*=0.1; e0222*=0.1;
            e1111*=0.1; e1112*=0.1; e1122*=0.1;
            e1212*=0.1; e1222*=0.1;
            e2222*=0.1;
            //obliczenie bledow iloczynow elementow tensora odkształceń
            double dE0000=0, dE0001=0, dE0002=0, dE0011=0, dE0012=0, dE0022=0,
                   dE0101=0, dE0102=0, dE0111=0, dE0112=0, dE0122=0,
                   dE0202=0, dE0211=0, dE0212=0, dE0222=0,
                   dE1111=0, dE1112=0, dE1122=0,
                   dE1212=0, dE1222=0,
                   dE2222=0;
            for (int i=0;i<10;i++) {double epsilon=e0000-e0000Set[i]; dE0000+=epsilon*epsilon;} dE0000=getAvErrorFromSumEps(dE0000,90.0); //10*9 (n(n-1))
            for (int i=0;i<10;i++) {double epsilon=e0001-e0001Set[i]; dE0001+=epsilon*epsilon;} dE0001=getAvErrorFromSumEps(dE0001,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0002-e0002Set[i]; dE0002+=epsilon*epsilon;} dE0002=getAvErrorFromSumEps(dE0002,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0011-e0011Set[i]; dE0011+=epsilon*epsilon;} dE0011=getAvErrorFromSumEps(dE0011,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0012-e0012Set[i]; dE0012+=epsilon*epsilon;} dE0012=getAvErrorFromSumEps(dE0012,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0022-e0022Set[i]; dE0022+=epsilon*epsilon;} dE0022=getAvErrorFromSumEps(dE0022,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0101-e0101Set[i]; dE0101+=epsilon*epsilon;} dE0101=getAvErrorFromSumEps(dE0101,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0102-e0102Set[i]; dE0102+=epsilon*epsilon;} dE0102=getAvErrorFromSumEps(dE0102,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0111-e0111Set[i]; dE0111+=epsilon*epsilon;} dE0111=getAvErrorFromSumEps(dE0111,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0112-e0112Set[i]; dE0112+=epsilon*epsilon;} dE0112=getAvErrorFromSumEps(dE0112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0122-e0122Set[i]; dE0122+=epsilon*epsilon;} dE0122=getAvErrorFromSumEps(dE0122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0202-e0202Set[i]; dE0202+=epsilon*epsilon;} dE0202=getAvErrorFromSumEps(dE0202,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0211-e0211Set[i]; dE0211+=epsilon*epsilon;} dE0211=getAvErrorFromSumEps(dE0211,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0212-e0212Set[i]; dE0212+=epsilon*epsilon;} dE0212=getAvErrorFromSumEps(dE0212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0222-e0222Set[i]; dE0222+=epsilon*epsilon;} dE0222=getAvErrorFromSumEps(dE0222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1111-e1111Set[i]; dE1111+=epsilon*epsilon;} dE1111=getAvErrorFromSumEps(dE1111,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1112-e1112Set[i]; dE1112+=epsilon*epsilon;} dE1112=getAvErrorFromSumEps(dE1112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1122-e1122Set[i]; dE1122+=epsilon*epsilon;} dE1122=getAvErrorFromSumEps(dE1122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1212-e1212Set[i]; dE1212+=epsilon*epsilon;} dE1212=getAvErrorFromSumEps(dE1212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1222-e1222Set[i]; dE1222+=epsilon*epsilon;} dE1222=getAvErrorFromSumEps(dE1222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e2222-e2222Set[i]; dE2222+=epsilon*epsilon;} dE2222=getAvErrorFromSumEps(dE2222,90.0);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... "); fflush(stdout);
            //Eijkl - bezwymiarowe [strain], volume - [sigma^3], T*=kT/eps [eps - jednostka energii (w domyśle =1)]
            double s0000=e0000*avVolume/TReduced, dS0000=(fabs(e0000*dAvVolume)+fabs(dE0000*avVolume))/TReduced,
                   s0001=e0001*avVolume/TReduced, dS0001=(fabs(e0001*dAvVolume)+fabs(dE0001*avVolume))/TReduced,
                   s0002=e0002*avVolume/TReduced, dS0002=(fabs(e0002*dAvVolume)+fabs(dE0002*avVolume))/TReduced,
                   s0011=e0011*avVolume/TReduced, dS0011=(fabs(e0011*dAvVolume)+fabs(dE0011*avVolume))/TReduced,
                   s0012=e0012*avVolume/TReduced, dS0012=(fabs(e0012*dAvVolume)+fabs(dE0012*avVolume))/TReduced,
                   s0022=e0022*avVolume/TReduced, dS0022=(fabs(e0022*dAvVolume)+fabs(dE0022*avVolume))/TReduced,
                   s0101=e0101*avVolume/TReduced, dS0101=(fabs(e0101*dAvVolume)+fabs(dE0101*avVolume))/TReduced,
                   s0102=e0102*avVolume/TReduced, dS0102=(fabs(e0102*dAvVolume)+fabs(dE0102*avVolume))/TReduced,
                   s0111=e0111*avVolume/TReduced, dS0111=(fabs(e0111*dAvVolume)+fabs(dE0111*avVolume))/TReduced,
                   s0112=e0112*avVolume/TReduced, dS0112=(fabs(e0112*dAvVolume)+fabs(dE0112*avVolume))/TReduced,
                   s0122=e0122*avVolume/TReduced, dS0122=(fabs(e0122*dAvVolume)+fabs(dE0122*avVolume))/TReduced,
                   s0202=e0202*avVolume/TReduced, dS0202=(fabs(e0202*dAvVolume)+fabs(dE0202*avVolume))/TReduced,
                   s0211=e0211*avVolume/TReduced, dS0211=(fabs(e0211*dAvVolume)+fabs(dE0211*avVolume))/TReduced,
                   s0212=e0212*avVolume/TReduced, dS0212=(fabs(e0212*dAvVolume)+fabs(dE0212*avVolume))/TReduced,
                   s0222=e0222*avVolume/TReduced, dS0222=(fabs(e0222*dAvVolume)+fabs(dE0222*avVolume))/TReduced,
                   s1111=e1111*avVolume/TReduced, dS1111=(fabs(e1111*dAvVolume)+fabs(dE1111*avVolume))/TReduced,
                   s1112=e1112*avVolume/TReduced, dS1112=(fabs(e1112*dAvVolume)+fabs(dE1112*avVolume))/TReduced,
                   s1122=e1122*avVolume/TReduced, dS1122=(fabs(e1122*dAvVolume)+fabs(dE1122*avVolume))/TReduced,
                   s1212=e1212*avVolume/TReduced, dS1212=(fabs(e1212*dAvVolume)+fabs(dE1212*avVolume))/TReduced,
                   s1222=e1222*avVolume/TReduced, dS1222=(fabs(e1222*dAvVolume)+fabs(dE1222*avVolume))/TReduced,
                   s2222=e2222*avVolume/TReduced, dS2222=(fabs(e2222*dAvVolume)+fabs(dE2222*avVolume))/TReduced,
                   //Voigt notation
                   S11=s0000, dS11=dS0000, S12=s0011, dS12=dS0011, S13=s0022, dS13=dS0022, S14=2*s0012, dS14=2*dS0012, S15=2*s0002, dS15=2*dS0002, S16=2*s0001, dS16=2*dS0001,
                   S22=s1111, dS22=dS1111, S23=s1122, dS23=dS1122, S24=2*s1112, dS24=2*dS1112, S25=2*s0211, dS25=2*dS0211, S26=2*s0111, dS26=2*dS0111,
                   S33=s2222, dS33=dS2222, S34=2*s1222, dS34=2*dS1222, S35=2*s0222, dS35=2*dS0222, S36=2*s0122, dS36=2*dS0122,
                   S44=4*s1212, dS44=4*dS1212, S45=4*s0212, dS45=4*dS0212, S46=4*s0112, dS46=4*dS0112,
                   S55=4*s0202, dS55=4*dS0202, S56=4*s0102, dS56=4*dS0102,
                   S66=4*s0101, dS66=4*dS0101,
                   //cubic average elements
                   S11c=(S11+S22+S33)/3.0, dS11c=(dS11+dS22+dS33)/3.0,
                   S12c=(S12+S13+S23)/3.0, dS12c=(dS12+dS13+dS23)/3.0,
                   S44c=(S44+S55+S66)/3.0, dS44c=(dS44+dS55+dS66)/3.0,
                   //isotropic elastic moduli
                   B=1/(3*S11c+6*S12c), dB=(fabs(dS11c)+2*fabs(dS12c))/(3*pow(S11c+2*S12c,2)),
                   my1=1/(2*S11c-2*S12c), dMy1=(fabs(dS11c)+fabs(dS12c))/(2*pow(S11c-S12c,2)),
                   my2=1/S44c, dMy2=fabs(dS44c/S44c/S44c),
                   avMy=(my1+my2)/2.0, dAvMy=(dMy1+dMy2)/2.0,
                   E=(9*B*avMy)/(3*B+avMy), dE=(9*(3*fabs(B*B*dAvMy)+fabs(dB*avMy*avMy)))/pow(3*B+avMy,2),
                   //Poisson ratio in key directions: [100], [111], [110][1m10], [110][001]
                   nu_100_all=-S12c/S11c, dNu_100_all=fabs(dS12c/S11c)+fabs((dS11c*S12c)/S11c/S11c),
                   nu_111_all=(-2*S11c-4*S12c+S44c)/(2*(S11c+2*S12c+S44c)), dNu_111_all=(3*(fabs(dS44c*(S11c+2*S12c))+(fabs(dS11c)+2*fabs(dS12c))*fabs(S44c)))/(2*pow(S11c+2*S12c+S44c,2)),
                   nu_110_1m10=1-(4*(S11c+S12c))/(2*(S11c+S12c)+S44c), dNu_110_1m10=(4*fabs(dS44c)*fabs(S11c+S12c)+4*(fabs(dS11c)+fabs(dS12c))*fabs(S44c))/pow(2*(S11c+S12c)+S44c,2),
                   nu_110_001=-((4*S12c)/(2*(S11c+S12c)+S44c)), dNu_110_001=(4*(2*fabs(dS11c)*fabs(S12c)+fabs(dS44c*S12c)+fabs(dS12c*(2*S11c+S44c))))/pow(2*(S11c+S12c)+S44c,2),
                   //Bij for stability conditions
                   B11=(S11c+S12c)/(S11c*S11c+S11c*S12c-2*S12c*S12c),   //>0
                   B44=1/S44c,                                          //>0
                   B12dB11=-S12c/(S11c+S12c);                           //>=-1/2 && <=1
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla 1 czastki
            /*printf("Creation of 1-particle orientation file for Origin... "); fflush(stdout);
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
            printf("done\n");*/

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla wszystkich czastek
            /*printf("Creation of ALL-particle orientation file for Origin... "); fflush(stdout);
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
            printf("done\n");*/

            //analiza konfiguracji przejsciowych (dPhi z konfiguracji na konfiguracje)
            /*double avAbsDPhi=0;
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
            }*/

            long timeEndMath=time(0);




/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... "); fflush(stdout);
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) fileConfigurations = fopen(bufferConfigurations,"w");
            fileOrientationsResults = fopen(bufferOrientationsResults,"w");
            fileAllOrientationsResults = fopen(allOrientationsResultsFileName,"w");

            fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avRho,dAvRho,avPacFrac,dAvPacFrac,avVolume,dAvVolume,avBoxMatrix[0],dAvBoxMatrix[0],avBoxMatrix[1],dAvBoxMatrix[1],avBoxMatrix[2],dAvBoxMatrix[2],avBoxMatrix[3],dAvBoxMatrix[3],avBoxMatrix[4],dAvBoxMatrix[4],avBoxMatrix[5],dAvBoxMatrix[5],B,dB,avMy,dAvMy,my1,dMy1,my2,dMy2,E,dE,nu_100_all,dNu_100_all,nu_111_all,dNu_111_all,nu_110_1m10,dNu_110_1m10,nu_110_001,dNu_110_001,S11c,dS11c,S12c,dS12c,S44c,dS44c,S11,dS11,S12,dS12,S13,dS13,S14,dS14,S15,dS15,S16,dS16,S22,dS22,S23,dS23,S24,dS24,S25,dS25,S26,dS26,S33,dS33,S34,dS34,S35,dS35,S36,dS36,S44,dS44,S45,dS45,S46,dS46,S55,dS55,S56,dS56,S66,dS66,B11,B44,B12dB11);
            fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",pressureReduced,avPacFrac,B,dB,avMy,dAvMy,E,dE,nu_100_all,dNu_100_all,nu_111_all,dNu_111_all,nu_110_1m10,dNu_110_1m10,nu_110_001,dNu_110_001,S11c,dS11c,S12c,dS12c,S44c,dS44c);

            if (!onlyMath[0]) {
                rho=N2/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2],deltaR,deltaPhi,deltaV,deltaVND);
                for (int i=0;i<activeN;i++) fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.17E,%.17E,%.17E,%.17E,%.17E,%.17E],",particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].phi[0],particles[i].phi[1],particles[i].phi[2],particles[i].halfSigma,particles[i].atoms[0]->diameter,particles[i].atoms[1]->diameter);
                //for (int i=0;i<activeN2;i++) fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E],",atoms[i].r[0],atoms[i].r[1],atoms[i].r[2]);  //aby tego uzyc, nalezy zakomentowac cialo metody 'computeAtomsPositions()' i uruchomic przy 0 cyklach rozruchowych i pomiarowych
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Parallelepiped[{0,0,0},{{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E}}]},{Opacity[0.2],Green,Sphere[{%.12E,%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][0],boxMatrix[1][1],boxMatrix[1][2],boxMatrix[2][0],boxMatrix[2][1],boxMatrix[2][2],particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].r[2],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[2][2]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E, boxMatrix[0][2]=boxMatrix[2][0]=%.12E, boxMatrix[1][2]=boxMatrix[2][1]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
            }

            //zakomentowane funkcje orientacyjne
            /*for (int i=0;i<ODFLength;i++) {
                fprintf(fileOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_1P[i]);
                fprintf(fileAllOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_AllP[i]);
            }*/

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) fclose(fileConfigurations);
            fclose(fileOrientationsResults); fclose(fileAllOrientationsResults);
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        getNextArgument(arg,true);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        oldTotalInteractionEnergy=totalInteractionEnergy;
        args[4]=0;
        loadedConfiguration=0;
        generatorStartPoint=0;
    }
    printf("\nTime for equilibrations: %ldsec, time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n"); fflush(stdout);
}
