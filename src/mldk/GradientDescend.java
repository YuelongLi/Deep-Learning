package mldk;

import java.util.Arrays;

/**
 * A Multi-Variable Optimization program that utilizes Yuelong's Gradient
 * Descend. It has great performance in most ranges, and costs minimal
 * iterations. The program demonstrates efficiency, accuracy, and flexibility.
 * <p>
 * Yuelong's gradient descend is an optimization method that uses dynamic
 * definition for dx and enhanced Newton's method.
 * </p>
 * 
 * @author Yuelong Li
 *
 */
public class GradientDescend {
	int size = 1;
	int iteratedTimes = 0;
	// the most time gradients can be performed
	int maxITS = 100;
	double lowestValue;
	/**
	 * dynamic definition of zero that changes as the coordinate approaches minimum
	 */
	private double ZERO = 0.000001;
	private double ZERO2 = 0.0001;
	final double Zero1;
	final double PSIX = 0.06;
	/**
	 * @param cds
	 *            current coordinates
	 */
	double[] cds;
	/**
	 * @param fps
	 *            Derivatives array
	 */
	double[] fps = new double[2];;
	double[] gradient;

	/**
	 * indicates whether the result is convergent and reliable
	 */
	public boolean convergent;

	/**
	 * Searches for the value closest to 0 when set to find zero mode, other wise it
	 * searches for the minimum
	 */
	boolean findZero;
	
	FunctionInterface fInterface;

	protected double function(double... coord) {
		return fInterface.function(coord);
	};
	
	public interface FunctionInterface{
		public double function(double... coords);
	}

	public GradientDescend(int maxITS, FunctionInterface... fInterface) {
		this.maxITS = maxITS;
		this.Zero1 = 0.000001;
		this.fInterface = fInterface[0];
	}

	public GradientDescend(int maxITS, double convergeRange,FunctionInterface... fInterface) {
		this.maxITS = maxITS;
		this.Zero1 = convergeRange;
		this.fInterface = fInterface[0];
	}

	/**
	 * finds the lowest point on the graph
	 * 
	 * @return whether the function is convergent
	 */
	public boolean find(double[] orgCoordinates) {
		this.cds = orgCoordinates.clone();
		findZero = false;
		gradient = new double[cds.length];
		for (int i = 0; i < maxITS; i++) {
			boolean converge = true;
			getGradient();
			Console.println("Gradients" + Arrays.toString(gradient));
			for (int a = 0; a < cds.length; a++) {
				double cd = cds[a];
				if (!Double.isNaN(gradient[a]))
					cd -= gradient[a];
				cds[a] = cd;
				if (checkValid()) {
					Console.println(a);
					while (checkValid()) {
						cd += 0.01 * Math.signum(gradient[a]);
						cds[a] = cd;
					}
					cds[a] = cd;
				}
				if (Math.abs(gradient[a]) > Zero1)
					converge = false;
			}

			Console.println("Iteration" + i);
			Console.println("cds" + Arrays.toString(cds));
			iteratedTimes++;
			if (converge)
				break;
		}
		return convergent;
	}

	public boolean findZero(double[] orgCoordinates) {
		this.cds = orgCoordinates.clone();
		findZero = true;
		gradient = new double[cds.length];
		for (int i = 0; i < maxITS; i++) {
			boolean converge = true;
			getGradient();
			Console.println("Gradients" + Arrays.toString(gradient));
			for (int a = 0; a < cds.length; a++) {
				double cd = cds[a];
				if (!Double.isNaN(gradient[a]))
					cd -= gradient[a];
				cds[a] = cd;
				if (checkValid()) {
					Console.println(a);
					while (checkValid()) {
						cd += 0.01 * Math.signum(gradient[a]);
						cds[a] = cd;
					}
					cd += 0.2 * Math.signum(gradient[a]);
					cds[a] = cd;
				}
				if (Math.abs(gradient[a]) > Zero1)
					converge = false;
			}

			Console.println("Iteration" + i);
			Console.println("cds" + Arrays.toString(cds));
			iteratedTimes++;
			if (converge)
				break;
		}
		return convergent;
	}

	protected void setZero(double numericSample) {
		int exponent = (int) Math.round(0.11536 * Math.log(100.58 + numericSample) - 6.8068);
		this.ZERO = Math.pow(10, exponent);
	}

	protected void setZero2(double numericSample) {
		int exponent = (int) Math.round(0.11536 * Math.log(100.58 + numericSample) - 6.8068);
		this.ZERO2 = Math.pow(10, exponent);
	}

	protected boolean checkValid() {
		return Double.isNaN(function(cds));
	}

	/**
	 * 
	 * @return the curved gradient at coordinate cds
	 * @param cds
	 *            coordinates of the current point
	 */
	private void getGradient() {
		convergent = true;
		setZero(function(cds));
		for (int i = 0; i < cds.length; i++) {
			fps(i);
			double curvedFDP = ((Math.abs(fps[1]) < 0.5) ? 0.5 : Math.abs(fps[1])) * ((findZero) ? Math.signum(fps[1]) : 1);
			Console.println("curvedFDP" + curvedFDP);
			try {
				gradient[i] = fps[0] / curvedFDP;
			} catch (ArithmeticException e) {
				Console.println("Can not optimize linear equation");
			}
			if ((int) fps[0] != 0)
				convergent = false;
		}
	}

	double ps;
	double p2s;

	/**
	 * finds the f prime f double prime and f triple prime of the assigned variable
	 * 
	 * @return big decimal value of first to the third derivatives
	 * @param variableIndex
	 *            specifies the variable of which the partial dirivative is computed
	 * @param coordinate
	 *            the coordinate of the variable specified
	 */
	private void fps(int variableIndex) {
		if (findZero) {
			ps = function(cds);
			p2s = getDeriv(cds, variableIndex);

			fps[0] = ps;
			fps[1] = p2s;
		} else {
			ps = getDeriv(cds, variableIndex);
			setZero2(ps);
			p2s = getDbDeriv(cds, variableIndex);
			fps[0] = ps;
			fps[1] = p2s;
		}
	}

	private double getDeriv(double[] cds, int variableIndex) {
		double[] cds1 = new double[cds.length];
		for (int i = 0; i < cds.length; i++) {
			cds1[i] = cds[i];
		}
		cds1[variableIndex] -= ZERO;
		return (function(cds) - function(cds1)) / ZERO;
	}

	private double getDbDeriv(double[] cds, int variableIndex) {
		double[] cds1 = new double[cds.length];
		for (int i = 0; i < cds.length; i++) {
			cds1[i] = cds[i];
		}
		cds1[variableIndex] -= ZERO2;
		return (getDeriv(cds, variableIndex) - getDeriv(cds1, variableIndex)) / ZERO2;
	}

	/**
	 * 
	 * @return the value of the function at the lowest point
	 */
	public double[] getResult() {
		double[] result = new double[cds.length + 1];
		for (int i = 0; i < cds.length; i++) {
			result[i] = cds[i];
		}
		result[cds.length] = function(cds);
		return result;
	}

	public String toString() {
		return "Descentd Result: " + Arrays.toString(getResult()) + ", Convergent: " + convergent + ", Iterated Times: "
				+ iteratedTimes + ", Gradient Accuracy: " + Zero1;

	}
}
/*
 * © Copyright 2016 Cannot be used without authorization
 */
