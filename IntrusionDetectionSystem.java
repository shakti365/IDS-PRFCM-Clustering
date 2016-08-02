/*
 * Novel Intrusion Detection System
 * The following project is divided into 2 phases and 5 major parts
 * Training Phase
 *  1. Intialisation of Cluster Centers
 *  2. Computation of nearest neighbours of trainset
 *  3. Penalization Reward based Fuzzy C means clustering
 *  4. Validity calculations
 * Testing Phase
 *  1. Modified weighted KNN classification
 *  2. Dempster-Shafer Theory
 * Input:
 *  1. KDD trainset
 *  2. KDD testset
 *  3. K nearest neighbours obtained from.........
 * Output:
 *  1. Confusion matrix for clustered trainset
 *  2. Confusion matrix for testset
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Scanner;

/**
 *
 * @author shakti
 */
public class IntrusionDetectionSystem {

    /**
     * @param args the command line arguments
     */
    static Double t[][] = new Double[125974][42];           //train set
    static Double c[][] = new Double[2][41];                //cluster center
    static Double m[][] = new Double[125974][2];            //membership matrix
    static Double m1[][] = new Double[125974][2];           //old membership matrix for iteration
    static Double v[][] = new Double[22545][42];            //test set
    static int truecounter = 0;
    static int falsecounter = 0;
    static int totiteration = 0;
    static int truepositive = 0;
    static int truenegative = 0;
    static int falsepositive = 0;
    static int falsenegative = 0;
    static int forj1 = 0;
    static int forj1still1 = 0;
    static int forj1true = 0;
    static int forj1false = 0;
    static int forj1notmatch = 0;
    static double rangeofweight;
    public static void main(String[] args) {
        // TODO code application logic here
        double w = 0.5d, max = 0.0d;
        int rows = 125974;      //number of rows in train dataset
        int columns = 41;       //number of columns in train dataset
        //train[][] stores the returned string array from getContentArray() 
        String train[][] = getContentArray("kddtrain_2class_normalized.csv", rows, columns + 1);
        //double t[][] array stores the converted train[][] array from getConvertedArray()
        t = getConvertedArray(train, rows, columns + 1);
        //for garbage collection
        for (int i = 0; i < rows; i++) {
            train[i] = null;
        }
        System.gc();
        int rowtest = 22545;    //number of rows in test dataset
        int coltest = 41;       //number of columns in test dataset
        //test[][] stores the returned string array from getContentArray() 
        String test[][] = getContentArray("KDDTest+_normalized_2.csv", rowtest, coltest + 1);
        //double v[][] array stores the converted test[][] array from getConvertedArray()
        v = getConvertedArray(test, rowtest, coltest + 1);
        //double c[][] array stores random cluster centers from getRandomInitialClusterCenter()
        c = getRandomInitialClusterCenter(t, rows, columns);
        //double c[][] array stores updated cluster centers by KNN from getKNNUpdatesClusterCenter()
        //c = getKNNUpdatedClusterCenter(c, rows, columns);
        //int p[][] array stores the returned neighbours from getNeighbours()
//        String position[][] = getContentArray("positionfinal_norm.csv", rows, 2);
        int pos[][] = getNeighbours(rows, 3);
//        position = null;
//        System.gc();
        //double m[][] array stores the initialized membership
        m = initMembership(t, c, rows, columns);
        int ctr = 0;
        do {
            System.out.println("\nITERATION No.: " + ctr++);

            //a1 & a2 values of alpha 1 and alpha 2
            double a1 = getAlpha(m, rows, 1);
            double a2 = getAlpha(m, rows, 2);

           // m = getPenaltyReward(m, pos, rows, 2);
            
//            System.out.println("\n====CLUSTER CENTER UPDATION====");            
            c = getUpdatedClusterCenter(m, t, rows, columns);

            for (int i = 1; i < rows; i++) {
                for (int j = 0; j < 2; j++) {
                    m1[i][j] = m[i][j];
                }
            }

           // m = initMembership(t, c, rows, columns);
            double maxi = 0.0d, mini = m[1][1];
            m = getUpdatedMembership(t, c, pos, m, rows, columns, a1, a2, w, maxi, mini);
            
            for (int i = 1; i < rows; i++) {
                for (int j = 0; j < 2; j++) {
                    m1[i][j] = Math.abs(m[i][j] - m1[i][j]);
                }

            }
            for (int i = 1; i < rows; i++) {
                for (int j = 0; j < 2; j++) {
                    if (m1[i][j] > max);
                    max = m1[i][j];
                }
            }
            
        }while(max > 0.0001);
        
        
        int output;
        for(int i=1;i<rows;i++){
            if(m[i][0]>m[i][1])
            {
                output = 1;
            }
            else{
                output = 2;
            }
            if((output == 1)&&(t[i][columns] == 1)){
                truepositive++;
                truecounter++;
                totiteration++;
            } else if((output == 1)&&(t[i][columns] == 2)){
                falsepositive++;
                falsecounter++;
                totiteration++;
            } else if((output == 2)&&(t[i][columns] == 2)){
                truenegative++;
                truecounter++;
                totiteration++;
            } else if((output == 2)&&(t[i][columns] == 1)){
                falsenegative++;
                falsecounter++;
                totiteration++;
            }
                
        }
        System.out.println(truepositive+"   "+falsepositive+"  "+truenegative+"  "+falsenegative);
        System.out.println(truecounter+"    "+falsecounter+"    "+totiteration);
        
        
        truecounter = 0;
        falsecounter = 0;
        totiteration = 0;
        truepositive = 0;
        truenegative = 0;
        falsepositive = 0;
        falsenegative = 0;
        
        setValidity(rows, 200);
        Double valid[] = getValidity(rows, 1);
        
        for(int whole=1; whole<rowtest; whole++){
            
            Double d[] = getDistanceConnection(whole, rows, columns);
            
            double min = d[1], max1 = 0.0d;
            int minpos=0;
            for(int i = 1; i < rows; i++){
                if (d[i] > max1) {
                max1 = d[i];
                }
                if (d[i] < min) {
                min = d[i];
                minpos = i;
                }
            }
            int j = 0;
            for (int i = 1; i < rows; i++) {
                if (d[i] == min) {
                    j++;
                }
            }
            
            if(j == 1){
                    forj1++;
                    Double newrule[][] =  getBestRuleSetFor1(d, whole, j, min, max1, valid, rows, columns, minpos);
                    if(newrule[0][0] != 2.0d){
                        output = Combination(newrule);
                        showResult(output, whole,1);
                    }
//                if(m[minpos][0]>m[minpos][1])
//                    output = 1;
//                else
//                    output = 2;
//                showResult(output, whole);
                
            } else {
                Double newrule[][] = getBestRuleSet(d, whole, j, min, max1, valid, rows, columns);
                output = Combination(newrule);
                showResult(output, whole,0);
            }
            System.out.println("range of weight: "+rangeofweight);
        }
        System.out.println("TOTAL TRUE: "+truecounter);
        System.out.println("TOTAL FALSE: "+falsecounter);
        System.out.println("TRUE positive: "+truepositive);
        System.out.println("TRUE negative: "+truenegative);
        System.out.println("false postive: "+falsepositive);
        System.out.println("false negative: "+falsenegative);
        System.out.println("TOTAL j=1: "+forj1);
        System.out.println("TOTAL j=1 still j=1: "+forj1still1);
        System.out.println("TOTAL j=1 TRUE: "+forj1true);
        System.out.println("TOTAL j=1 FALSE: "+forj1false);
        System.out.println("number of time minpos is not selected: "+forj1notmatch);
    }

    /**
     * getContentArray() takes a comma separated value(.csv) file and returns
     * its content as an array Arguments: fileName: name of the file rows:
     * number of rows columns: number of columns Return: content: a
     * String[rows][columns] array of the .csv file
     */
    private static String[][] getContentArray(String fileName, int rows, int columns) {
        String content[][] = new String[rows][columns];
        try {
            Scanner s = new Scanner(new File(fileName));
            int i = 0;
            while (s.hasNext()) {

                String temp = s.nextLine();
                content[i++] = temp.split(",", columns + 1);
            }
        } catch (FileNotFoundException e) {
            System.err.println("file not found");
        }
        return content;
    }



    /**
     * getConvertedArray() takes a comma separated value(.csv) file and returns
     * its content as an array Arguments: StringArray: string array to be
     * converted to double rows: number of rows columns: number of columns
     * Return: DoubleArray: a Double[rows][columns] converted array of the
     * StringArray
     */
    private static Double[][] getConvertedArray(String StringArray[][], int rows, int columns) {
        Double DoubleArray[][] = new Double[rows][columns];
        for (int i = 1; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                DoubleArray[i][j] = Double.parseDouble(StringArray[i][j]);
            }
        }
        return DoubleArray;
    }
    private static int[][] getConvertedArrayToInt(String StringArray[][], int rows, int columns) {
        int IntegerArray[][] = new int[rows][columns];
        for (int i = 1; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                IntegerArray[i][j] = Integer.parseInt(StringArray[i][j]);
            }
        }
        return IntegerArray;
    }

    /**
     * getInitialClusterCenter() takes the train set and returns randomly
     * generated cluster centers Arguments: TrainArray: Double array containing
     * train set to find cluster from rows: number of rows columns: number of
     * columns Return: ClusterCenter: Double[2][columns] array containing
     * initial randomly selected clusters
     */
    private static Double[][] getRandomInitialClusterCenter(Double TrainArray[][], int rows, int columns) {
        Double ClusterCenter[][] = new Double[2][columns];
        for (int j = 0; j < columns; j++) {
            ClusterCenter[0][j] = TrainArray[rows - 2][j];
            ClusterCenter[1][j] = TrainArray[rows - 1][j];
        }
        return ClusterCenter;
    }

    /**
     * getKNNUpdatedClusterCenter() returns the updated cluster centers by KNN
     * method Arguments: c[][]: cluster center rows: number of rows columns:
     * number of columns Return: c[][]: updated cluster center
     */
    private static Double[][] getKNNUpdatedClusterCenter(Double c[][], int rows, int columns) {
        int blahblah = 0;
        do {
            int ctr1 = 0, ctr2 = 0;   //counter of nearest neighbours of clusters 1 & 2
            Double d1[][] = new Double[rows][3];
            //calculating distance of each connection to the intial cluser center
            for (int i = 1; i < rows; i++) {
                d1[i][0] = getDistance(t, i, c, 0, columns);
                d1[i][1] = getDistance(t, i, c, 1, columns);
                if (d1[i][0] < d1[i][1]) //decision of the connection which cluster it belongs 1 or 2
                {
                    d1[i][2] = 1.0;
                } else {
                    d1[i][2] = 2.0;
                }

            }
            for (int i = 1; i < rows; i++) {
                if (d1[i][2] == 1.0) {
                    ctr1++;
                    for (int j = 0; j < columns; j++) {
                        c[0][j] += t[i][j];
                    }
                }
                if (d1[i][2] == 2.0) {
                    ctr2++;
                    for (int j = 0; j < columns; j++) {
                        c[1][j] += t[i][j];
                    }
                }

            }
            for (int j = 0; j < columns; j++) {   //assigning new clusters
                c[0][j] = (double) c[0][j] / ctr1;
                c[1][j] = (double) c[1][j] / ctr2;
            }
            blahblah++;
        } while (blahblah != 3);    // cluster center intialisation completed
        return c;
    }

    /**
     * getNeighbours() `returns the position of K-nearest neighbours Arguments:
     * rows: number of rows k: number of neighbours Return: pos: pos[rows][2]
     * array containing K nearest neighbour positions
     */
    private static int[][] getNeighbours(int rows, int k) {
        String position[][] = getContentArray("position_10_norm.csv", rows, k);
        int pos[][] = new int[rows][k];
        for (int i = 1; i < rows; i++) {
            for (int j = 0; j < k; j++) {
                double x = Double.parseDouble(position[i][j]);
                pos[i][j] = (int) x;
            }
            position[i] = null;
        }
        System.gc();
        return pos;
    }

    /**
     * getDistance() returns the distance between two vectors Parameters: t[][]:
     * train set i: row of train set c[][]: cluster center array j: row of
     * cluster center dimension: dimension of the vector Return: c[][]
     * initialized cluster centers
     */
    private static double getDistance(Double t[][], int i, Double c[][], int j, int dimension) {
        double dis1 = 0.0d;
        for (int k = 0; k < dimension; k++) {
            dis1 += Math.pow((t[i][k] - c[j][k]), 2.0d);
        }
        return (Math.sqrt(dis1));
    }

    /**
     * initMembership() returns the initial membership Parameters: t[][]: train
     * set c[][]: cluster center array j: row columns: columns Return: m[][]
     * initialized membership
     */
    private static Double[][] initMembership(Double t[][], Double c[][], int rows, int columns) {
        Double m[][] = new Double[rows][2];
        for (int i = 1; i < rows; i++) {
            double denom1, denom2;  //variables for membership calculation
            denom1 = getDistance(t, i, c, 0, columns);
            denom2 = getDistance(t, i, c, 1, columns);
            m[i][0] = 1.0d / (1.0d + Math.pow((denom1 / denom2), 2.0d));
            m[i][1] = 1.0d / (Math.pow((denom2 / denom1), 2.0d) + 1.0d);
//            System.out.println(m[i][0] + "\t" + m[i][1]);     //display membership
        }
        return m;
    }

    /**
     * getAlpha() returns the value of alpha 1 and alpha 2 Parameters: m[][]:
     * membership array rows: row index: alpha index(1 or 2) Return: alpha:
     * calculated alpha formula
     */
    private static double getAlpha(Double m[][], int rows, int index) {
        double q1 = 0.0d, q2 = 0.0d;
        for (int i = 1; i < rows; i++) {
            
                q1 += Math.pow(m[i][0], 2);
            
                q2 += Math.pow(m[i][1], 2);
        }
        if (index == 1) {
            return (q1 / (q1 + q2));
        } else if (index == 2) {
            return (q2 / (q1 + q2));
        }
        return 0;
    }

    /**
     * getUpdatedClusterCenter() returns the updated value of cluster center
     * using FCM Parameters: m[][]: membership array t[][]: training set rows:
     * row columns: columns Return: c[][]: cluster center
     */
    private static Double[][] getUpdatedClusterCenter(Double m[][], Double t[][], int rows, int columns) {
        double x, y;
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < columns; j++) {
                double num = 0.0d;
                double denom = 0.0d;
                for (int i = 1; i < rows; i++) {
                    x = Math.pow((m[i][k]), 2.0d);
                    y = t[i][j];
                    num = num + (x * y);
                    denom = denom + x;
                }
                c[k][j] = num / denom;
            }
        }
        return c;
    }

    
    /**
     * getUpdatedMembership()    returns the updated value of membership matrix using PRFCM 
     * Parameters:
     
     *  t[][]:          training set
     *  c[][]:          cluster center
     *  pos[][]:        nearest neighbour positions
     *  m[][]:          membership array
     *  rows:           row 
     *  columns:        columns
     *  a1:             alpha1
     *  a2:             alpha2
     * Return:
     *  m[][]:          membership matrix
     */
    private static Double[][] getUpdatedMembership(Double t[][], Double c[][], int pos[][], Double m[][], int rows, int columns, double a1, double a2, double w, double maximum, double minimum){
        for (int i = 1; i <rows; i++) {
                double den1 = 0.0, den2 = 0.0, denom1, denom2;  //variables for membership calculation
                for (int j = 0; j < columns; j++) {
                    den1 += Math.pow((t[i][j] - c[0][j]), 2.0d);
                    den2 += Math.pow((t[i][j] - c[1][j]), 2.0d);
                }
                Double reward[][] = getRewardPenalty(w, m, pos, i, a1, a2, 3);
                m[i][0] = (Math.pow((den1 + reward[0][0] - reward[0][1]), -1)) / (Math.pow((den1 + reward[0][0] - reward[0][1]), -1) + Math.pow((den2 + reward[1][1] - reward[1][0]), -1));
                m[i][1] = (Math.pow((den2 + reward[1][1] - reward[1][0]), -1)) / (Math.pow((den1 + reward[0][0] - reward[0][1]), -1) + Math.pow((den2 + reward[1][1] - reward[1][0]), -1));
//                double eps = 0.000000001d;
//                if(m[i][0] < 0.0d)
//                    m[i][0] = eps;
//                if(m[i][1] < 0.0d)
//                    m[i][1] = eps;
//                if(m[i][0] > 1.0d)
//                    m[i][0] = 1.0d - eps; 
//                if(m[i][1] > 1.0d)
//                    m[i][1] = 1.0d - eps;
//                
                if(m[i][0] < minimum)
                    minimum = m[i][0];
                if(m[i][1] < minimum)
                    minimum = m[i][1];
                if(m[i][0] > maximum)
                    maximum = m[i][0];
                if(m[i][1] > maximum)
                    maximum = m[i][1];
//                    m[i][0] = Math.abs(m[i][0]) / (Math.abs(m[i][0]) + Math.abs(m[i][1]) );
//                    m[i][1] = Math.abs(m[i][1]) / (Math.abs(m[i][0]) + Math.abs(m[i][1]) );
                
//                m[i][0] = (Math.pow((den1 - reward[0][0] + reward[1][0]), -1)) / (Math.pow((den1 - reward[0][0] + reward[1][0]), -1) + Math.pow((den2 - reward[0][1] + reward[1][1]), -1));
//                m[i][1] = (Math.pow((den2 + reward[1][1] - reward[0][1]), -1)) / (Math.pow((den1 + reward[1][0] - reward[0][0]), -1) + Math.pow((den2 + reward[1][1] - reward[0][1]), -1));

                
//                m[i][0] = (Math.pow((den1 + reward[0][0]), -1)) / (Math.pow((den1 + reward[0][0]), -1) + Math.pow((den2 + reward[0][1]), -1));
//                m[i][1] = (Math.pow((den2 + reward[1][1]), -1)) / (Math.pow((den1 + reward[1][0]), -1) + Math.pow((den2 + reward[1][1]), -1));
             
        
        }  
        m = normalize(m, maximum, minimum, rows);
        return m;
    }
    private static Double[][] normalize(Double m[][], double max, double min, int rows){
        for(int i = 1; i < rows; i++){
            if(m[i][0] == min)
                m[i][0] = 0.00000001d;
            else if(m[i][0] == max)
                m[i][0] = 1.0d - 0.00000001d;
            else
                m[i][0] = (m[i][0] - min) / (max - min);
            m[i][1] = 1.0d - m[i][0];
        }
        return m;
    }
    /**
     * getUpdatedMembership()    returns the updated value of membership matrix using PRFCM 
     * Parameters:
     
     *  t[][]:          training set
     *  c[][]:          cluster center
     *  pos[][]:        nearest neighbour positions
     *  m[][]:          membership array
     *  rows:           row 
     *  columns:        columns
     *  a1:             alpha1
     *  a2:             alpha2
     * Return:
     *  m[][]:          membership matrix
     */
    private static Double[][] getRewardPenalty(double w, Double m[][], int pos[][], int t, double a1, double a2, int k){
    
        Double p[][] = new Double[k][2];
        Double reward[][] = new Double[2][2];
        reward[0][0] = 0.0d;
        reward[0][1] = 0.0d;
        reward[1][0] = 0.0d;
        reward[1][1] = 0.0d;
        for(int i = 0; i < k; i++){
            p[i][0] = m[pos[t][i]][0];
            p[i][1] = m[pos[t][i]][1];
        }
        for(int i = 0; i < k; i++){
            reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
            reward[0][1] += (Math.pow(-1, 2) * w * (1.0d - p[i][0]) * Math.log(1.0d - a1));
            reward[1][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
            reward[1][0] += (Math.pow(-1, 2) * w * (1.0d - p[i][1]) * Math.log(1.0d - a2));
//            if(m[t][0]>m[t][1]){
//                if(p[i][0] > p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * (1.0d - p[i][0]) * Math.log(a2));
//                    reward[1][1] += (Math.pow(-1, 2) * w * (1.0d - p[i][1]) * Math.log(a1));
//                }
//                else if(p[i][0] < p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a2));
//                    reward[1][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a1));
//                }
//            }
//            else if(m[t][0]<m[t][1]){
//                if(p[i][0] < p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a2));
//                    reward[1][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a1));
//                }
//                else if(p[i][0] > p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a2));
//                    reward[1][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a1));
//                }
//            }
            
//            if(m[t][0]>m[t][1]){
//                if(p[i][0] > p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * (1.0d - p[i][0]) * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 2) * w * (1.0d - p[i][1]) * Math.log(a2));
//                }
//                else if(p[i][0] < p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * (1.0d - p[i][0]) * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 2) * w * (1.0d - p[i][1]) * Math.log(a2));
//                }
//            }
//            else if(m[t][0]<m[t][1]){
//                if(p[i][0] < p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * (1.0d - p[i][0]) * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 2) * w * (1.0d - p[i][1]) * Math.log(a2));
//                }
//                else if(p[i][0] > p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * (1.0d - p[i][0]) * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 2) * w * (1.0d - p[i][1]) * Math.log(a2));
//                }
//            }
            
//            if(m[t][0]>m[t][1]){
//                if(p[i][0] > p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 1) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 1) * w * p[i][0] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                }
//                else if(p[i][0] < p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 1) * w * p[i][1] * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 1) * w * p[i][1] * Math.log(a2));
//                }
//            }
//            else if(m[t][0]<m[t][1]){
//                if(p[i][0] < p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 1) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 1) * w * p[i][0] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 2) * w * p[i][1] * Math.log(a2));
//                }
//                else if(p[i][0] > p[i][1]){
//                    reward[0][0] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a1));
//                    reward[0][1] += (Math.pow(-1, 2) * w * p[i][0] * Math.log(a2));
//                    reward[1][0] += (Math.pow(-1, 1) * w * p[i][1] * Math.log(a1));
//                    reward[1][1] += (Math.pow(-1, 1) * w * p[i][1] * Math.log(a2));
//                }
//            }
        }
        return reward;
    }
    private static void setValidity(int rows, int k){
        Double valid[] = new Double[rows];
        for(int i = 1; i < 10000; i++){
            String position[] = getPos("positionfinal1_norm.csv", k, i);
            int s = 0;
            for(int j = 0; j < k-1; j++){
                double x= Double.parseDouble(position[j]);
                int y = (int)x;
                if(((m[y][0]>m[y][1])&&(m[i][0]>m[i][1]))||((m[y][0]<m[y][1])&&(m[i][0]<m[i][1]))){
                    s += 1;
                } else {
                    s += 0;
                }
            }
            valid[i] = (double)s/k;
            System.out.println(i);
        }
        writeFile("validity.csv", valid, 10000);
    }
    public static void writeFile(String filename, Double valid[], int rows){
        PrintWriter p = null;
        try{
            File file = new File(filename);
            p = new PrintWriter(new FileWriter(file, true));
            for(int i = 1; i < rows; i++){
                p.println(valid[i]);
            }
        } catch(Exception e){
            System.err.println("no file");
        } finally{
            p.close();
        }
        
    }
    private static String[] getPos(String fileName, int columns, int n) { //reads and returns the content of a .csv file in an array
        String content[] = new String[columns];
        try {
            Scanner s = new Scanner(new File(fileName));
           
            for(int i = 0; i < n; i++)
                s.nextLine();
            String temp = s.nextLine();
            content = temp.split(",", columns+1);
            
        } catch (Exception e) {
            System.err.println("file not found");
        }
        return content;
    }
    private static Double[] getValidity(int rows, int columns){
        String validity[][] = getContentArray("validity.csv", rows, columns);
        Double valid[] = new Double[rows];
        for(int i = 1; i < rows; i++){
            valid[i] = Double.parseDouble(validity[i][0]);
        }
        return valid;
    }
    
    private static Double[] getDistanceConnection(int whole, int rows, int columns){
        Double d[] = new Double[rows];
        for (int i = 1; i < rows; i++) {
            double sqd = 0.0d;
            for (int j = 0; j < columns; j++) {
                sqd += Math.pow((t[i][j] - v[whole][j]), 2.0d);     //connection used here
            }
            d[i] = Math.sqrt(sqd);
        }
        return d;
    }
    
    
    private static Double[][] getBestRuleSet(Double d[], int whole, int j, double min, double max1, Double valid[], int rows, int columns){
        Double b[][] = new Double[j][columns];
        int pos[] = new int[j];
        int q = 0;
        for(int i = 1; i < rows; i++){
            if (d[i] == min) {
                pos[q] = i;
                for (int k = 0; k < columns; k++) {
                    b[q][k] = t[i][k];
                }   
                q++;
            }
        }
        
        Double weight[] = new Double[j];
        for(int i=0; i<j; i++){
            weight[i] = valid[pos[i]] * ((double)(max1-d[pos[i]])/(max1-min));
        }
        
        Double con[][] = new Double[rows][columns];
        for(int i=1; i<rows; i++){
            double temp =0.0d;
            for(j=0; j<2; j++){
                temp = temp + m[i][j];
            }
            double p = 1.0d/temp;
            for(j=0; j<2; j++){
                con[i][j] = p * m[i][j];
            }
        }
        
        //RULESET MATRIX
        Double r[][] = new Double[rows][columns];
        int count=0;
        for(int i=0; i<pos.length; i++){
            for(j=0; j<2; j++){
                double rule = weight[i] * con[pos[i]][j];
                if(rule != 0.0){
                    r[i][j] = rule;
                    count++;
                }
            }
        }
        
        //NEW RULESET
        Double newrule[][] = new Double[pos.length][3];
        for(int i=0;i<pos.length;i++){
            double sum=0.0d;
            for(j=0;j<2;j++){
                
                if(r[i][j]!=0.0d){
                    newrule[i][j] = r[i][j];
                    sum += newrule[i][j];
                }
            }
            newrule[i][2] = 1.0d - sum;
        }
        
        return newrule;
    }
    
    private static int Combination(Double newrule[][]){
        //COMBINATION
        int output=0;
        double kinv = 1.0d/(1.0d-((newrule[0][0]*newrule[1][1])+(newrule[0][1]*newrule[1][0])));
        double mc1 = ((newrule[0][0]*newrule[1][0])+(newrule[0][0]*newrule[1][2])+(newrule[0][2]*newrule[1][0]))*kinv;
        double mc2 = ((newrule[0][1]*newrule[1][1])+(newrule[0][1]*newrule[1][2])+(newrule[0][2]*newrule[1][1]))*kinv;
        double mun = (newrule[0][2]*newrule[1][2])*kinv;
        double belc1 =mc1;
        double belc2 =mc2;
        double plc1 = 1.0d - belc1;
        double plc2 = 1.0d - belc2;
        double unc1 = belc1 - plc1;
        double unc2 = belc2 - plc2;
        double bpc1 = mc1 + (mun/2);
        double bpc2 = mc2 + (mun/2);
        if(bpc1>bpc2)
            output = 1;
        else
            output = 2;
        return output;
    }
    
    private static void showResult(int output, int whole, int flag){
        if(flag == 1){
            System.out.println("out: "+output+"\t test: "+v[whole][41]);
            if(output == v[whole][41]){
                System.out.println("true: "+(++truecounter)+" jtrue: "+(forj1true++));
            }
            else
                System.out.println("false: "+(++falsecounter)+" jfalse: "+(forj1false++));
            System.out.println("total: "+(++totiteration));
        }
        else{
            System.out.println("out: "+output+"\t test: "+v[whole][41]);
            if(output == v[whole][41]){
                System.out.println("true: "+(++truecounter));
            }
            else
                System.out.println("false: "+(++falsecounter));
            System.out.println("total: "+(++totiteration));
        }
            
        if(output == 1 && v[whole][41] == 1)
            System.out.println("truepositive: "+(++truepositive));
        else if(output == 2 && v[whole][41] == 1)
            System.out.println("falsepositive: "+(++falsepositive));
        else if(output == 2 && v[whole][41] == 2)
            System.out.println("truenegative: "+(++truenegative));
        else if(output == 1 && v[whole][41] == 2)
            System.out.println("falsenegative: "+(++falsenegative));
        
    }
    private static Double[][] getBestRuleSetFor1(Double d[], int whole, int j, double min, double max1, Double valid[], int rows, int columns, int minpos){
        double max = 0.0d;
        Double weight[] = new Double[rows];
        int count=0;
        for(int i = 1; i < rows; i++){
//            if(i == minpos)
//                continue;
            weight[i] = valid[i] * ((double)(max1-d[i])/(max1-min));
            if(weight[i]>max){
                max = weight[i];
            }
        }
        for(int i =1 ; i <rows;i++){
            if(weight[i] == max)
                count++;
        }
        
        //System.out.println(count);    
        int pos[] = new int[count];
        Double b[][] = new Double[count][columns];
        j = 0;
        
        for(int i = 1; i < rows; i++){
            if(max == weight[i]){
                pos[j] = i;
                for(int l = 0; l < columns; l++){
                    b[j][l] = t[i][l]; 
                }
                j++;
            }
        }
        int flag = 0;
       
        for(int i = 0; i < j; i++){
            if(pos[i] == minpos){
               
                flag = 1;
            }
        }
        rangeofweight = Math.abs(weight[minpos] - weight[pos[0]]);
        if(flag == 0)
            forj1notmatch++;
        int output=0;
        if(count == 1){
            forj1still1++;
            if(m[pos[0]][0] > m[pos[0]][1])
                    output = 1;
            else 
                    output = 2;
                showResult(output, whole,1);
                b[0][0] = 2.0d;
                 return b;
        }
        else{
        
        Double con[][] = new Double[rows][columns];
        for(int i=1; i<rows; i++){
            double temp =0.0d;
            for(j=0; j<2; j++){
                temp = temp + m[i][j];
            }
            double p = 1.0d/temp;
            for(j=0; j<2; j++){
                con[i][j] = p * m[i][j];
            }
        }
        
        //RULESET MATRIX
        Double r[][] = new Double[rows][columns];
        count=0;
        for(int i=0; i<pos.length; i++){
            for(j=0; j<2; j++){
                double rule = weight[pos[i]] * con[pos[i]][j];
                if(rule != 0.0){
                    r[i][j] = rule;
                    count++;
                }
            }
        }
        
        //NEW RULESET
        Double newrule[][] = new Double[pos.length][3];
        for(int i=0;i<pos.length;i++){
            double sum=0.0d;
            for(j=0;j<2;j++){
                if(r[i][j]!=0.0d){
                    newrule[i][j] = r[i][j];
                    sum += newrule[i][j];
                }
            }
            newrule[i][2] = 1.0d - sum;
        }
        
        return newrule;
    }
       
    }
    
    private static Double[][] getPenaltyReward(Double m[][], int pos[][], int rows, int k){
        double x=0.0d,y=0.0d;
        for(int i = 1; i < rows; i++){
            for(int j = 0; j < k; j++){
                if(m[i][0] > m[i][1]){
                    if(m[pos[i][j]][0] > m[pos[i][j]][1]){
                        x += m[pos[i][j]][0];
                        y -= m[pos[i][j]][1];
                    }
                    else if(m[pos[i][j]][0] < m[pos[i][j]][1]){
                        x -= m[pos[i][j]][0];
                        y += m[pos[i][j]][1];
                    }
                }
                else if(m[i][0] < m[i][1]){
                    if(m[pos[i][j]][0] < m[pos[i][j]][1]){
                        x += m[pos[i][j]][0];
                        y -= m[pos[i][j]][1];
                    }
                    else if(m[pos[i][j]][0] > m[pos[i][j]][1]){
                        x -= m[pos[i][j]][0];
                        y += m[pos[i][j]][1];
                    }
                }
            }
            m[i][0] = x/(x+y);
            m[i][1] = y/(x+y);
        }
        return m;
    }
    
}
