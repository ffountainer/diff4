// packages for the graph
package org.example;
import javax.swing.*;
import java.util.ArrayList;


import org.jzy3d.chart.Chart;
import org.jzy3d.chart.ChartLauncher;
import org.jzy3d.chart.factories.AWTChartComponentFactory;
import org.jzy3d.chart.factories.IChartComponentFactory;
import org.jzy3d.chart.factories.SwingChartComponentFactory;
import org.jzy3d.maths.Coord3d;
import org.jzy3d.plot3d.primitives.Scatter;
import org.jzy3d.plot3d.rendering.canvas.ICanvas;
import org.jzy3d.plot3d.rendering.canvas.Quality;


public class Main {
    public static void main(String[] args) {
        solve("Task", 0, 1, 0, 0);
    }

    public static void solve(String task, double x_0, double y_0, double z_0, double t_0_) {
        ArrayList<ArrayList<Triple>> secondOrderAnswers = new ArrayList<>();
        ArrayList<ArrayList<Double>> time = new ArrayList<>();
        double[] deltaT = {0.001, 0.0005};
        for (int delta_index = 0; delta_index < deltaT.length; delta_index++) {
            secondOrderAnswers.add(new ArrayList<>());

            double x = x_0;
            double y = y_0;
            double z = z_0;
            double t = t_0_;

            ArrayList<Double> timeCur = new ArrayList<>();

            while (t <= 100.01) {
                Quadruple quad = secondOrderRungeKuttaIteration(y, x, z, t, delta_index, 1, deltaT);
                t = quad.firstElement;
                x = quad.secondElement;
                y = quad.thirdElement;
                z = quad.fourthElement;
                Triple triple = new Triple(x, y, z);
                timeCur.add(t);
                secondOrderAnswers.get(delta_index).add(triple);

            }

            time.add(timeCur);

        }


        ArrayList<ArrayList<Double>> secondOrderErrors = new ArrayList<>();


        for (int i = 1; i < secondOrderAnswers.size(); i++) {
            secondOrderErrors.add(new ArrayList<>());

            int step2 = (int) ((int)Math.round(100/deltaT[i]));
            int step1 = (int) ((int)Math.round(100/deltaT[i - 1]));

            int ratio2 = step2 / 100;
            int ratio1 = step1 / 100;

            int k = ratio2;
            int j = ratio1;


            for (int counter = 0; counter < 100; counter++) {

                double elementx1 = secondOrderAnswers.get(i - 1).get(j).getFirstElement();
                double elementx2 = secondOrderAnswers.get(i).get(k).getFirstElement();


                double elementy1 = secondOrderAnswers.get(i - 1).get(j).getSecondElement();
                double elementy2 = secondOrderAnswers.get(i).get(k).getSecondElement();

                double elementz1 = secondOrderAnswers.get(i - 1).get(j).getThirdElement();
                double elementz2 = secondOrderAnswers.get(i).get(k).getThirdElement();


                secondOrderErrors.get(i - 1).add(Math.abs(elementx2 - elementx1) + Math.abs(elementy2 - elementy1) + Math.abs(elementz2 - elementz1));

                j += ratio1;
                k += ratio2;
            }
        }

        int ind = 0;
        int num = -1;
        for (ArrayList<Double> errors : secondOrderErrors) {
            double absError = errors.get(0);
            if (absError < 1e-8) {
                num = ind;
            }
            for (double error : errors) {
                if (error > absError && error <  1e-8) {
                    num = ind;
                }
                absError = Math.max(absError, error);
            }

            ind += 1;
        }

        int[] N = new int[deltaT.length];
        for (int i = 0; i< deltaT.length; i++) {
            N[i] = (int) ((int)Math.round(100/deltaT[i]));
        }

        System.out.println();
        System.out.println("Second order");

        if (num == -1) {
            num = 0;
            System.out.println("The max error is more than 10^-8");
        }

        System.out.println("Local errors for 2nd order Runge-Kutta:");
        System.out.println("N = " + N[num]);
        for (double error : secondOrderErrors.get(num)) {
            System.out.print(error + " ");
        }

        System.out.println();

        System.out.println();


        ArrayList<ArrayList<Triple>> fourthOrderAnswers = new ArrayList<>();

        for (int delta_index = 0; delta_index < deltaT.length; delta_index++) {
            fourthOrderAnswers.add(new ArrayList<>());

            double x = x_0;
            double y = y_0;
            double z = z_0;
            double t = t_0_;

            ArrayList<Double> timeCur = new ArrayList<>();

            while (t <= 100.01) {
                Quadruple quad = fourthOrderRungeKuttaIteration(y, x, z, t, delta_index, 1, deltaT);
                t = quad.firstElement;
                x = quad.secondElement;
                y = quad.thirdElement;
                z = quad.fourthElement;
                Triple triple = new Triple(x, y, z);
                timeCur.add(t);
                fourthOrderAnswers.get(delta_index).add(triple);

            }

            time.add(timeCur);

        }


        ArrayList<ArrayList<Double>> fourthOrderErrors = new ArrayList<>();


        for (int i = 1; i < fourthOrderAnswers.size(); i++) {
            fourthOrderErrors.add(new ArrayList<>());

            int step2 = (int) ((int)Math.round(100/deltaT[i]));
            int step1 = (int) ((int)Math.round(100/deltaT[i - 1]));

            int ratio2 = step2 / 100;
            int ratio1 = step1 / 100;

            int k = ratio2;
            int j = ratio1;


            for (int counter = 0; counter < 100; counter++) {

                double elementx1 = fourthOrderAnswers.get(i - 1).get(j).getFirstElement();
                double elementx2 = fourthOrderAnswers.get(i).get(k).getFirstElement();


                double elementy1 = fourthOrderAnswers.get(i - 1).get(j).getSecondElement();
                double elementy2 = fourthOrderAnswers.get(i).get(k).getSecondElement();

                double elementz1 = fourthOrderAnswers.get(i - 1).get(j).getThirdElement();
                double elementz2 = fourthOrderAnswers.get(i).get(k).getThirdElement();


                fourthOrderErrors.get(i - 1).add(Math.abs(elementx2 - elementx1) + Math.abs(elementy2 - elementy1) + Math.abs(elementz2 - elementz1));

                j += ratio1;
                k += ratio2;
            }
        }

        int ind4 = 0;
        int num4 = -1;
        for (ArrayList<Double> errors : fourthOrderErrors) {
            double absError = errors.get(0);
            if (absError < 1e-8) {
                num4 = ind4;
            }
            for (double error : errors) {
                if (error > absError && error <  1e-8) {
                    num4 = ind4;
                }
                absError = Math.max(absError, error);
            }

            ind += 1;
        }


        int[] N4 = new int[deltaT.length];
        for (int i = 0; i< deltaT.length; i++) {
            N4[i] = (int) ((int)Math.round(100/deltaT[i]));
        }


        System.out.println("Fourth order");

        if (num4 == -1) {
            num4 = 0;
            System.out.println("The max error is more than 10^-8");
        }

        System.out.println("Local errors for 4nd order Runge-Kutta:");
        System.out.println("N = " + N[num4]);
        for (double error : fourthOrderErrors.get(num4)) {
            System.out.print(error + " ");
        }


        System.out.println();



        double[] secondOrderDataY = new double[secondOrderAnswers.get(num).size()];
        for (int i = 0; i < secondOrderAnswers.get(num).size(); i++) {
            secondOrderDataY[i] = secondOrderAnswers.get(num).get(i).getSecondElement();
        }


        double[] secondOrderDataX = new double[secondOrderAnswers.get(num).size()];
        for (int i = 0; i < secondOrderAnswers.get(num).size(); i++) {
            secondOrderDataX[i] = secondOrderAnswers.get(num).get(i).getFirstElement();
        }

        double[] secondOrderDataZ = new double[secondOrderAnswers.get(num).size()];
        for (int i = 0; i < secondOrderAnswers.get(num).size(); i++) {
            secondOrderDataZ[i] = secondOrderAnswers.get(num).get(i).getThirdElement();
        }

        Chart chart1 = SwingChartComponentFactory.chart(Quality.Advanced, IChartComponentFactory.Toolkit.swing);

        Coord3d[] points1 = new Coord3d[secondOrderDataX.length];
        for (int i = 0; i < secondOrderDataX.length; i++) {
            points1[i] = new Coord3d(secondOrderDataX[i], secondOrderDataY[i], secondOrderDataZ[i]);
        }

        Scatter scatter1 = new Scatter(points1);

        chart1.getScene().add(scatter1);

        ICanvas canvas1 = chart1.getCanvas();

        JFrame frame1 = new JFrame("3D Plot Example");
        frame1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        frame1.getContentPane().add((java.awt.Component) canvas1);

        frame1.setSize(800, 600);
        frame1.setVisible(true);


        double[] fourthOrderDataY = new double[fourthOrderAnswers.get(num4).size()];
        for (int i = 0; i < fourthOrderAnswers.get(num4).size(); i++) {
            fourthOrderDataY[i] = fourthOrderAnswers.get(num4).get(i).getSecondElement();
        }

        double[] fourthOrderDataX = new double[fourthOrderAnswers.get(num4).size()];
        for (int i = 0; i < fourthOrderAnswers.get(num4).size(); i++) {
            fourthOrderDataX[i] = fourthOrderAnswers.get(num4).get(i).getFirstElement();
        }

        double[] fourthOrderDataZ = new double[fourthOrderAnswers.get(num4).size()];
        for (int i = 0; i < fourthOrderAnswers.get(num4).size(); i++) {
            fourthOrderDataZ[i] = fourthOrderAnswers.get(num4).get(i).getThirdElement();
        }

        Chart chart = SwingChartComponentFactory.chart(Quality.Advanced, IChartComponentFactory.Toolkit.swing);

        Coord3d[] points = new Coord3d[fourthOrderDataX.length];
        for (int i = 0; i < fourthOrderDataX.length; i++) {
            points[i] = new Coord3d(fourthOrderDataX[i], fourthOrderDataY[i], fourthOrderDataZ[i]);
        }

        Scatter scatter = new Scatter(points);

        chart.getScene().add(scatter);

        ICanvas canvas = chart.getCanvas();

        JFrame frame = new JFrame("3D Plot Example");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        frame.getContentPane().add((java.awt.Component) canvas); // Correct type cast

        frame.setSize(800, 600);
        frame.setVisible(true);
    }

    public static Quadruple secondOrderRungeKuttaIteration(double y_t, double x_t, double z_t, double t, int deltaTInd, int direction, double[] deltaT) {

        double h = deltaT[deltaTInd];

        double dy_1 = 1.0;
        double dy_2 = 1.0;

        double dx_1 = 1.0;
        double dx_2 = 1.0;

        double dz_1 = 1.0;
        double dz_2 = 1.0;

        dx_1 = -3.0 * (x_t - y_t);
        dy_1 = 26.5 * x_t - y_t - x_t * z_t;
        dz_1 = x_t * y_t - z_t;

        dx_2 = -3.0 * ((x_t + (h / 2.0) * dx_1) - (y_t + (h / 2.0) * dy_1));
        dy_2 = 26.5 * (x_t + (h / 2.0) * dx_1) - (y_t + (h / 2.0) * dy_1) - (x_t + (h / 2.0) * dx_1) * (z_t + (h / 2.0) * dz_1);
        dz_2 = (x_t + (h / 2.0) * dx_1) * (y_t + (h / 2.0) * dy_1) - (z_t + (h / 2.0) * dz_1);

        x_t = x_t + (dx_1 + dx_2) * 0.5;
        y_t = y_t + (dy_1 + dy_2) * 0.5;
        z_t = z_t + (dz_1 + dz_2) * 0.5;

        t = t + h;

        return new Quadruple(t, x_t, y_t, z_t);

    }

    public static Quadruple fourthOrderRungeKuttaIteration(double y_t, double x_t, double z_t, double t, int deltaTInd, int direction, double[] deltaT) {
        deltaT[deltaTInd] = deltaT[deltaTInd] * direction;

        double dy_1 = 1.0;
        double dy_2 = 1.0;
        double dy_3 = 1.0;
        double dy_4 = 1.0;

        double dx_1 = 1.0;
        double dx_2 = 1.0;
        double dx_3 = 1.0;
        double dx_4 = 1.0;

        double dz_1 = 1.0;
        double dz_2 = 1.0;
        double dz_3 = 1.0;
        double dz_4 = 1.0;

        double h = deltaT[deltaTInd];

        dx_1 = -3 * (x_t - y_t);
        dy_1 = 26.5 * x_t - y_t - x_t * z_t;
        dz_1 = x_t * y_t - z_t;

        dx_2 = -3 * ((x_t + (h/2) * dx_1) - (y_t + (h/2) * dy_1));
        dy_2 = 26.5 * (x_t + (h/2) * dx_1) - (y_t + (h/2) * dy_1) - (x_t + (h/2) * dx_1) * (z_t + (h/2) * dz_1);
        dz_2 = (x_t + (h/2) * dx_1) * (y_t + (h/2) * dy_1) - (z_t + (h/2) * dz_1);

        dx_3 = -3 * ((x_t + (h/2) * dx_2) - (y_t + (h/2) * dy_2));
        dy_3 = 26.5 * (x_t + (h/2) * dx_2) - (y_t + (h/2) * dy_2) - (x_t + (h/2) * dx_2) * (z_t + (h/2) * dz_2);
        dz_3 = (x_t + (h/2) * dx_2) * (y_t + (h/2) * dy_2) - (z_t + (h/2) * dz_2);

        dx_4 = -3 * ((x_t + h * dx_3) - (y_t + h * dy_3));
        dy_4 = 26.5 * (x_t + h * dx_3) - (y_t + h * dy_3) - (x_t + h * dx_3)*(z_t + h * dz_3);
        dz_4 = (x_t + h * dx_3) * (y_t + h * dy_3) - (z_t + h * dz_3);

        x_t = x_t + (h/6) * (dx_1 + 2 * dx_2 + 2 * dx_3 + dx_4);
        y_t = y_t + (h/6) * (dy_1 + 2 * dy_2 + 2 * dy_3 + dy_4);
        z_t = z_t + (h/6) * (dz_1 + 2 * dz_2 + 2 * dz_3 + dz_4);

        t = t + h;


        return new Quadruple(t, x_t, y_t, z_t);

    }

    private static double roundTo10(double in) {
        int integer = (int) (in * Math.pow(10, 10));
        return (double) integer;
    }
}

class Quadruple{
    double firstElement;
    double secondElement;
    double thirdElement;
    double fourthElement;

    public Quadruple(double first, double second, double third, double fourth) {
        this.firstElement = first;
        this.secondElement = second;
        this.thirdElement = third;
        this.fourthElement = fourth;

    }

    public double getFirstElement() {
        return firstElement;
    }

    public double getSecondElement() {
        return secondElement;
    }

    public double getThirdElement() {
        return thirdElement;
    }

    public double getFourthElement() {
        return fourthElement;
    }
}

class Pair{
    double firstElement;
    double secondElement;

    public Pair(double first, double second) {
        this.firstElement = first;
        this.secondElement = second;

    }

    public double getFirstElement() {
        return firstElement;
    }

    public double getSecondElement() {
        return secondElement;
    }

}

class Triple{
    double firstElement;
    double secondElement;
    double thirdElement;

    public Triple(double first, double second, double third) {
        this.firstElement = first;
        this.secondElement = second;
        this.thirdElement = third;

    }

    public double getFirstElement() {
        return firstElement;
    }

    public double getSecondElement() {
        return secondElement;
    }

    public double getThirdElement() {
        return thirdElement;
    }
}
