import ilog.concert.*;
import ilog.cplex.*;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

public class Cplex {

    static double[][] distanceMatrix, probMatrix;
    static int noOfVillages;

    public static void partA() {
        try {
            IloCplex cplex = new IloCplex();

            //x[j]: 1 if village j is chosen; 0 otherwise
            IloNumVar[] x = cplex.boolVarArray(noOfVillages);

            IloNumVar[][] y = new IloNumVar[noOfVillages][];
            for (int i = 0; i < noOfVillages; i++)
                y[i] = cplex.boolVarArray(noOfVillages);

            IloNumVar T = cplex.numVar(0.0, Double.MAX_VALUE);

            cplex.addMinimize(T);


            //yij <= xj for all i,j
            for (int i = 0; i < noOfVillages; i++)
                for (int j = 0; j < noOfVillages; j++)
                    cplex.addLe(cplex.sum(y[i][j], cplex.prod(-1.0, x[j])), 0);

            //sum of xi's = 4
            cplex.addEq(cplex.sum(x), 4);


            //sum of yij over j = 1 for all i
            for (int i = 0; i < noOfVillages; i++)
                cplex.addEq(cplex.sum(y[i]), 1);


            //sum of dij*yij over j less then T
            for (int i = 0; i < noOfVillages; i++)
                cplex.addLe(cplex.sum(cplex.scalProd(distanceMatrix[i], y[i]), cplex.prod(-1.0, T)), 0);

             // solve the model and display the solution if one was found
             if ( cplex.solve() ) {
                 double[] xx = cplex.getValues(x);
                 double[][] yy = new double[noOfVillages][];

                 for(int i = 0; i < noOfVillages; i++)
                    yy[i] = cplex.getValues(y[i]);

                 System.out.println("\n\nPART A\n\"Solution status = " + cplex.getStatus() + "\nMax distance = " + cplex.getObjValue());
                 System.out.print("Villages chosen by santa: ");
                 for (int j = 0; j < noOfVillages; j++)
                    if(xx[j] == 1.0)
                        System.out.print( j + "  ");

                 System.out.println("\n-----------------------------------------\n");
            }
             cplex.end();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void partB(){
        try {
            IloCplex cplex = new IloCplex();

            //x[j]: 1 if village j is chosen; 0 otherwise
            IloNumVar[] x = cplex.boolVarArray(noOfVillages);

            IloNumVar[][] y = new IloNumVar[noOfVillages][];
            for(int i = 0; i < noOfVillages; i++)
                y[i] = cplex.boolVarArray(noOfVillages);

            IloNumVar T = cplex.numVar(0.0, Double.MAX_VALUE);

            cplex.addMinimize(T);


            //yij <= xj for all i,j
            for(int i=0; i < noOfVillages; i++)
                for (int j = 0; j < noOfVillages; j++)
                    cplex.addLe( cplex.sum( y[i][j], cplex.prod( -1.0, x[j]) ), 0 );

            //sum of xi's = 4
            cplex.addEq( cplex.sum(x), 4);


            //sum of yij over j = 1 for all i
            for(int i = 0; i < noOfVillages; i++)
                cplex.addEq(cplex.sum(y[i]), 1);


            //sum of dij*yij over j less then T
            for(int i = 0; i < noOfVillages; i++)
                cplex.addLe( cplex.sum( cplex.scalProd(distanceMatrix[i], y[i]), cplex.prod(-1.0, T) ), 0 );


            // pij>0.60 implies yij=0
            IloNumVar[][] a = new IloNumVar[noOfVillages][];
            for(int i = 0; i < noOfVillages; i++)
                a[i] = cplex.boolVarArray(noOfVillages);
            
            for(int i=0; i < noOfVillages; i++) {
                for (int j = 0; j < noOfVillages; j++) {
                    cplex.addLe( cplex.sum(  probMatrix[i][j] - 0.6,  cplex.prod(-1, a[i][j]) ), 0 );
                    cplex.addLe( cplex.sum( -1, cplex.sum( y[i][j], a[i][j]) ),0);
                }
            }



            if ( cplex.solve() ) {
                double[] xx = cplex.getValues(x);

                System.out.println("\n\nPART B\nSolution status =" + cplex.getStatus() + "Max Distance: " + cplex.getObjValue());
                System.out.print("Villages chosen by santa: ");
                for (int j = 0; j < noOfVillages; j++) {
                    if (xx[j] == 1.0)
                        System.out.print(j + "  ");

                }

                System.out.println("\n-----------------------------------------\n");

            }

            cplex.end();
        }
        catch (IloException e) { e.printStackTrace(); }
    }

    public static void partC(){
        try {
            IloCplex cplex = new IloCplex();

            //x[i][j]: 1 if santa travels from i directly to j
            IloNumVar[][] x = new IloNumVar[noOfVillages][];
            for(int i = 0; i < noOfVillages; i++)
                x[i] = cplex.boolVarArray(noOfVillages);

            IloNumVar T = cplex.numVar(0.0, Double.MAX_VALUE);

            //minimize total travel time
            cplex.addMinimize( T );

            IloNumExpr sum = cplex.sum( 0.0, cplex.scalProd(distanceMatrix[0], x[0]) );
            for(int i = 1; i < noOfVillages; i++){
                sum = cplex.sum( sum, cplex.scalProd(distanceMatrix[i], x[i]) );
            }
            sum = cplex.prod(0.025, sum);
            cplex.addEq(T, sum);

            // sum of xij over j's = 1
            for(int i = 0; i < noOfVillages; i++)
                cplex.addEq(cplex.sum(x[i]), 1);

            //sum of xij over i's = 1
            for(int j = 0; j < noOfVillages; j++) {
                sum = cplex.sum( 0, x[0][j] );
                for (int i = 1; i < noOfVillages; i++)
                    sum = cplex.sum( sum, x[i][j] );
                cplex.addEq( sum, 1);
            }

            for(int i = 0; i < noOfVillages; i++)
                cplex.addEq( x[i][i], 0);



            //t[j]>=t[i]+1 - (2*no_of_locs)*(1-x[i,j])
            //tj -ti -1 + 2*N >= + 2*Nxij

            //t: time at which village i is visited
            IloNumVar[] t = new IloNumVar[noOfVillages]; //cplex.numVarArray(noOfVillages, 1, noOfVillages-1);

            t[0] = cplex.numVar(0,0);

            for(int i=1; i < noOfVillages; i++)
                t[i] = cplex.numVar(1,noOfVillages-1);


            for(int i=1; i < noOfVillages; i++) {
                for (int j = 1; j < noOfVillages; j++) {
                    if(i != j)
                        cplex.addLe(    cplex.sum(1,cplex.diff(t[i], t[j]))    , cplex.prod(noOfVillages-1, cplex.diff(1, x[i][j])));
                }
            }


            // pij>0.60 implies xij=0
            IloNumVar[][] a = new IloNumVar[noOfVillages][];
            for(int i = 0; i < noOfVillages; i++)
                a[i] = cplex.boolVarArray(noOfVillages);


            for(int i=0; i < noOfVillages; i++) {
                for (int j = 0; j < noOfVillages; j++) {
                    cplex.addLe( cplex.sum(  probMatrix[i][j] - 0.6,  cplex.prod(-1, a[i][j]) ), 0 );
                    cplex.addLe( cplex.sum( -1, cplex.sum( x[i][j], a[i][j]) ),0);
                }
            }




            if ( cplex.solve() ) {
                double[][] xx = new double[noOfVillages][];

                for(int i = 0; i < noOfVillages; i++)
                    xx[i] = cplex.getValues(x[i]);

                //double[] tt = cplex.getValues(t);


                System.out.println("\n");
                System.out.println("PART C");

                cplex.output().println("Solution status = " + cplex.getStatus());
                cplex.output().println("Solution value = " + cplex.getObjValue());

                System.out.println("Total travel time: " + cplex.getValue(T));
                System.out.print("Path: ");

                int count = 0, i=0;
                System.out.print(i + " > ");
                while(count < noOfVillages){
                    for (int j = 0; j < noOfVillages; j++) {
                        if(xx[i][j] == 1) {
                            System.out.print(j + " > ");
                            i = j;
                            break;
                        }
                    }
                    count++;
                }


                System.out.println();
                System.out.println("-----------------------------------------");
                System.out.println();

            }

            cplex.end();
        }
        catch (IloException e) { e.printStackTrace(); }
    }

    public static void partD(){
        try {
            IloCplex cplex = new IloCplex();

            //x[i][j][k]: 1 if volunteer k travels from i to j
            IloNumVar[][][] x = new IloNumVar[noOfVillages][noOfVillages][];
            for(int i = 0; i < noOfVillages; i++)
                for(int j = 0; j < noOfVillages; j++)
                    x[i][j] = cplex.boolVarArray(noOfVillages);

            IloNumVar[] v = cplex.boolVarArray(noOfVillages);

            //# of volunteers
            cplex.addMinimize(cplex.sum(v));

            IloNumExpr sum1, sum2, sum3, sum4,sum5;
            //co1
            sum1 = cplex.prod(0, x[0][0][0]);
            for(int j = 1; j < noOfVillages; j++)
                cplex.sum(sum1, cplex.sum(x[0][j]));
            cplex.addEq(sum1, cplex.sum(v), "c1");

            //co2
            sum2 = cplex.prod(0, x[0][0][0]);
            for(int i = 1; i < noOfVillages; i++)
                cplex.sum(sum2, cplex.sum(x[i][0]));
            cplex.addEq(sum2, cplex.sum(v), "co2");


            //co3-i
            for(int i = 1; i < noOfVillages; i++) {
                sum3 = cplex.prod(0, x[0][0][0]);
                for (int j = 0; j < noOfVillages; j++) {
                    sum3 = cplex.sum(sum3, cplex.sum(x[i][j]));
                }
                cplex.addEq(sum3,1, "co3-"+i);
            }

            //co4-j
            for(int j = 1; j < noOfVillages; j++) {
                sum4 = cplex.prod(0, x[0][0][0]);
                for (int i = 0; i < noOfVillages; i++) {
                    cplex.sum(sum4, cplex.sum(x[i][j]));
                }
                cplex.addEq(sum4, 1,"co4-"+j);
            }

            //co5-i
            for(int i = 0; i < noOfVillages; i++)
                cplex.addEq( cplex.sum(x[i][i]), 0, "co5-"+i);



            for(int k = 0; k< noOfVillages; k++){
                sum5 = cplex.prod(0, x[0][0][0]);
                for(int i = 0; i < noOfVillages; i++){
                    for(int j = 0; j < noOfVillages; j++){
                        cplex.sum(sum5, cplex.prod(distanceMatrix[i][j], x[i][j][k], v[k]));
                    }
                }
                cplex.addLe(sum5, 400);
            }



            //t: time at which village i is visited
            IloNumVar[] t = new IloNumVar[noOfVillages]; //cplex.numVarArray(noOfVillages, 1, noOfVillages-1);

            t[0] = cplex.numVar(0,0);

            for(int i=1; i < noOfVillages; i++)
                t[i] = cplex.numVar(1,noOfVillages-1);


            for(int i=1; i < noOfVillages; i++) {
                for (int j = 1; j < noOfVillages; j++) {
                    if(i != j)
                        cplex.addLe(    cplex.sum(1,cplex.diff(t[i], t[j]))    , cplex.prod(noOfVillages-1, cplex.diff(1, cplex.sum(x[i][j]))));
                }
            }


            if ( cplex.solve() ) {
                double[][][] xx = new double[noOfVillages][][];

                for(int i = 0; i < noOfVillages; i++)
                    for(int j = 0; j < noOfVillages; j++)
                        xx[i][j] = cplex.getValues(x[i][j]);

                double[] vv = cplex.getValues(v);


                System.out.println("\n");
                System.out.println("PART C");

                cplex.output().println("Solution status = " + cplex.getStatus());
                cplex.output().println("Solution value = " + cplex.getObjValue());

                //System.out.println("# of volunteers: " + cplex.getValue(V));
                System.out.print("Path: ");


                for(int k = 0; k < noOfVillages && vv[k] == 1; k++){

                    int i  =0;
                    System.out.print(i + " > ");
                    for (int j = 0; j < noOfVillages; j++) {
                        if(xx[i][j][k] == 1) {
                            System.out.print(j + " > ");
                            i = j;
                            break;
                        }
                    }
                }



                System.out.println();
                System.out.println("\n-----------------------------------------\n");
                System.out.println();

            }

            cplex.end();
        }
        catch (IloException e) { e.printStackTrace(); }
    }


    public static void main(String[] args) throws IOException {

        //String fileName = "C:\\Users\\Tayyar Önal\\Desktop\\RABİA\\4.1\\IE400\\proje\\Cplex\\src\\data.xlsx";
        String fileName = args[0];

        File file = new File(fileName);
        FileInputStream fis = new FileInputStream(file);
        XSSFWorkbook wb = new XSSFWorkbook(fis);

        XSSFSheet sheetD = wb.getSheetAt(0);
        XSSFSheet sheetP = wb.getSheetAt(1);

        noOfVillages = sheetD.getPhysicalNumberOfRows();
        distanceMatrix = new double[noOfVillages][noOfVillages];
        probMatrix = new double[noOfVillages][noOfVillages];

        for(int i=0; i < noOfVillages; i++)
            for(int j=0; j < noOfVillages; j++) {
                distanceMatrix[i][j] = sheetD.getRow(i).getCell(j).getNumericCellValue();
                probMatrix[i][j] = sheetP.getRow(i).getCell(j).getNumericCellValue();
            }

        partA();
        partB();
        partC();
        partD();

    }

}
