
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class resultmerger {
    String N,gaps,G,mN,mS,mD,badanieStart,badanieEnd,newDirectoryName;
    
    String getDirectoryName(String badanieNumber) {
        return "2D_N-"+N+"_gaps-"+gaps+"_G-"+G+"_badanie-"+badanieNumber+"_mN-"+mN+"_mS-"+mS+"_mD-"+mD;
    }
    
    public resultmerger (String newDirectoryName, String N, String gaps, String G, String mN, String mS, String mD, String badanieStart, String badanieEnd, int liczbaDanych) {
        this.newDirectoryName=newDirectoryName;
        this.N=N;
        this.gaps=gaps;
        this.G=G;
        this.mN=mN;
        this.mS=mS;
        this.mD=mD;
        this.badanieStart=badanieStart;
        this.badanieEnd=badanieEnd;
           
        File resultsDirectory = new File(newDirectoryName); resultsDirectory.mkdir();
        resultsDirectory = new File(newDirectoryName+"/"+getDirectoryName(String.valueOf(badanieStart))); resultsDirectory.mkdir();
          
        for (int i=Integer.parseInt(badanieStart);i<=Integer.parseInt(badanieEnd);i++) {
            File actualDirectory = new File(getDirectoryName(String.valueOf(i)));
            for (int j=0;j<actualDirectory.list().length;j++) {
                if (actualDirectory.list()[j].endsWith("Results.txt")) {
                    try {
                        BufferedReader czytnik = new BufferedReader(new FileReader(actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]));
                        int skipChars = actualDirectory.list()[j].split("_")[0].length();
                        BufferedWriter zapis = new BufferedWriter(new FileWriter(resultsDirectory.getAbsolutePath()+"/j-none_"+actualDirectory.list()[j].substring(skipChars+1),true));
                        String buffer; while ((buffer=czytnik.readLine())!=null) {
                            try {
                                String[] dataS = buffer.split("\t");
                                for (int k=0;k<=liczbaDanych;k++) Double.parseDouble(dataS[k]);
                                zapis.write(buffer); zapis.newLine();
                            } catch (Exception e3) {System.out.println("Bad line found - fixed in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]);}
                        }
                        zapis.close(); czytnik.close();
                    } catch (Exception e2) {System.err.print("Error in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]+". "); e2.printStackTrace();}
                } else if (actualDirectory.list()[j].endsWith("transient.txt")) {
                    try {
                        BufferedReader czytnik = new BufferedReader(new FileReader(actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]));
                        int skipChars = actualDirectory.list()[j].split("_")[0].length();
                        BufferedWriter zapis = new BufferedWriter(new FileWriter(resultsDirectory.getAbsolutePath()+"/j-none_"+actualDirectory.list()[j].substring(skipChars+1),true));
                        zapis.write("newCfgFile"); zapis.newLine();
                        String buffer; while ((buffer=czytnik.readLine())!=null) {
                            zapis.write(buffer); zapis.newLine();
                        }
                        zapis.close(); czytnik.close();
                    } catch (Exception e2) {System.err.print("Error in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]+". "); e2.printStackTrace();}
                } else if (actualDirectory.list()[j].endsWith("allOnt.txt") && actualDirectory.list()[j].contains("Orientations")) {
                    try {//sprawdza pierwszy wiersz, gdy zawiera on wiecej ujemnych katow (50%), to laczenie nastepuje z odwroceniem katow. W efekcie, wyzszy pik zawsze jest dodatni.
                        BufferedReader czytnik = new BufferedReader(new FileReader(actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]));
                        String[] negativeAngleCheck = czytnik.readLine().split(",");
                        int negativeAngleCounter = 0;
                        for (int k=0;k<negativeAngleCheck.length;k++) if (Double.parseDouble(negativeAngleCheck[k])<0) negativeAngleCounter++;
                        czytnik.close();
                        
                        czytnik = new BufferedReader(new FileReader(actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]));
                        int skipChars = actualDirectory.list()[j].split("_")[0].length();
                        BufferedWriter zapis = new BufferedWriter(new FileWriter(resultsDirectory.getAbsolutePath()+"/j-none_"+actualDirectory.list()[j].substring(skipChars+1),true));
                        String buffer; 
                        double negDivAllAngle = (double)negativeAngleCounter/(double)negativeAngleCheck.length;
                        if (negDivAllAngle>0.5) {
                            System.out.println("Domination of negative orientations - converting in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]);
                            while ((buffer=czytnik.readLine())!=null) {
                                try {
                                    String[] bufferValueTable = buffer.split(",");
                                    String bufferString=String.valueOf(Double.parseDouble(bufferValueTable[0])*(-1));
                                    for (int k=1;k<bufferValueTable.length;k++) {
                                        zapis.write(bufferString);
                                        bufferString=String.valueOf(Double.parseDouble(bufferValueTable[k])*(-1));
                                        zapis.write(",");
                                    }
                                    zapis.write(bufferString);
                                    zapis.newLine();
                                } catch (Exception e3) {System.err.print("Error in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]+". "); e3.printStackTrace();}
                            }
                        } else 
                            while ((buffer=czytnik.readLine())!=null) {
                                zapis.write(buffer); zapis.newLine();
                            }
                        zapis.close(); czytnik.close();
                    } catch (Exception e2) {System.err.print("Error in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]+". "); e2.printStackTrace();}
                } else if (actualDirectory.list()[j].contains("Orientations")) {
                    try {
                        BufferedReader czytnik = new BufferedReader(new FileReader(actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]));
                        int skipChars = actualDirectory.list()[j].split("_")[0].length();
                        File zapisFile = new File(resultsDirectory.getAbsolutePath()+"/j-none_"+actualDirectory.list()[j].substring(skipChars+1)); 
                        boolean nieIstnieje = true; if (zapisFile.exists()) nieIstnieje = false;
                        BufferedWriter zapis = new BufferedWriter(new FileWriter(zapisFile,true));
                        if (nieIstnieje) zapis.append("{"); else zapis.append(",");
                        try {
                            String buffer = czytnik.readLine(); buffer=buffer.substring(1,buffer.length()-1);
                            zapis.append(buffer); 
                        } catch (Exception e3) {}
                        zapis.close(); czytnik.close();
                    } catch (Exception e2) {System.err.print("Error in: "+actualDirectory.getAbsolutePath()+"/"+actualDirectory.list()[j]+". "); e2.printStackTrace();}
                }
            }
        }
        for (int i=0;i<resultsDirectory.list().length;i++) {
            try {
                if (resultsDirectory.list()[i].contains("Orientations") && !resultsDirectory.list()[i].endsWith("allOnt.txt")) {
                    BufferedWriter zapis = new BufferedWriter(new FileWriter(resultsDirectory.getAbsolutePath()+"/"+resultsDirectory.list()[i],true));
                    zapis.append("}"); zapis.close();
                }
            } catch (Exception e2) {System.err.print("Error in: "+resultsDirectory.getAbsolutePath()+"/"+resultsDirectory.list()[i]+". "); e2.printStackTrace();}
        }
    }

    public static void main(String[] args) {
        new resultmerger(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],Integer.parseInt(args[9]));
    }
}
