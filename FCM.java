import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Scanner;

class FCM{
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
    public static void main(String[] args) {
    	double w = 0.0d, max = 0.0d;
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
//        int pos[][] = getNeighbours(rows, 3);
//        position = null;
//        System.gc();
        //double m[][] array stores the initialized membership
        m = initMembership(t, c, rows, columns);
        int ctr = 0;
        do {
            System.out.println("\nITERATION No.: " + ctr++);

            //a1 & a2 values of alpha 1 and alpha 2
            //double a1 = getAlpha(m, rows, 1);
            //double a2 = getAlpha(m, rows, 2);

           // m = getPenaltyReward(m, pos, rows, 2);
            
//            System.out.println("\n====CLUSTER CENTER UPDATION====");            
            c = getUpdatedClusterCenter(m, t, rows, columns);

            for (int i = 1; i < rows; i++) {
                for (int j = 0; j < 2; j++) {
                    m1[i][j] = m[i][j];
                }
            }

            m = initMembership(t, c, rows, columns);
            //double maxi = 0.0d, mini = m[1][1];
            //m = getUpdatedMembership(t, c, pos, m, rows, columns, a1, a2, w, maxi, mini);
            
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
                if(m[i][0] < minimum)
                    minimum = m[i][0];
                if(m[i][1] < minimum)
                    minimum = m[i][1];
                if(m[i][0] > maximum)
                    maximum = m[i][0];
                if(m[i][1] > maximum)
                    maximum = m[i][1];
             
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
        }
        return reward;
    }
}