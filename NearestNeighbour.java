/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Scanner;

/**
 *
 * @author shakti
 */
public class NearestNeighbour2 {
    public static void main(String[] args) {
        int rows = 125974;//number of rows from data set    dataset k w
        int columns = 41;//number of columns from dataset
        Double t[][] = new Double[rows][columns+1];    //train set 
        
        String train[][] = getContentArray("kddtrain_2class_normalized.csv", rows, columns+1);
        for (int i = 1; i < rows; i++) {
            for (int j = 0; j < columns+1; j++) {
                t[i][j] = Double.parseDouble(train[i][j]);
            }
        }
        
        for(int i = 1; i < rows; i++){
            System.out.println(i);
            Double d[][] = new Double[rows][2];
            double p[] = new double[rows];
            double dis1 = 0.0d;
            double min = 10000.0d;
            int m = 0;
            for(int j = 1; j < rows; j++){
                if(j != i){
                    double ds1 = 0.0d;
                    for(int k = 0; k < columns; k++){
                        ds1 += Math.pow((t[i][k] - t[j][k]), 2.0d);
                    }
                    d[m][0] = Math.sqrt(ds1);
                    d[m][1] = (double)j;
                    m++;
                }
            }

            for(int j = 0; j < 20; j++){
                for(int k = j+1; k < rows-2; k++){
                    if(d[j][0] > d[k][0]){
                        double x;
                        x = d[j][0];
                        d[j][0] = d[k][0];
                        d[k][0] = x;
                        x = d[j][1];
                        d[j][1] = d[k][1];
                        d[k][1] = x;
                    }
                }
            }
            int z=0;
            for(int j = 0; j < 20; j++){
                Double de = new Double(d[j][1]);
                int q = de.intValue();
                System.out.println(q);
                if (t[i][columns] == t[q][columns]){
                    System.out.println(t[i][columns] +" "+ t[q][columns]);
                
                p[z++] = d[j][1];
                }
            }
            writeFile("position_10_norm.csv", p, i);
        }
        
    }
    private static String[][] getContentArray(String fileName, int rows, int columns) { //reads and returns the content of a .csv file in an array
        String content[][] = new String[rows][columns];
        try {
            Scanner s = new Scanner(new File(fileName));
            int i = 0;
            while (s.hasNext()) {
                String temp = s.nextLine();
                content[i++] = temp.split(",");
            }
        } catch (Exception e) {
            System.err.println("file not found");
        }
        return content;
    }
    
    public static void writeFile(String filename, double x[], int j){
        PrintWriter p = null;
        try{
            File file = new File(filename);
            p = new PrintWriter(new FileWriter(file, true));
            for(int i = 1; i < 2; i++){
                if(j==i)
                    continue;
                p.print(x[i]+",");
            }
            p.println();
            
        } catch(Exception e){
            System.err.println("no file");
        } finally{
            p.close();
        }
        
    }
    
}
