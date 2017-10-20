package mldk;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

/**
 * The class integrates many operation regarding the console, 
 * such as reading from console, writing from console, etc.
 * @author Yuelong Li
 *
 */

public class Console {
	
	private static Scanner sc = new Scanner(System.in);
	
	public static void println(){
		System.out.println();
	}
	
	public static void println(Object a){
		if(a instanceof double[]) System.out.println(Arrays.toString((double[]) a));
		else System.out.println(a);
	};
	
	public static void print(Object a){
		if(a instanceof double[]) System.out.print(Arrays.toString((double[]) a));
		else System.out.print(a);
	}
	
	public static void println(int[] a){
		Console.println(Arrays.toString(a));
	}
	public static void print(int[] a){
		Console.print(Arrays.toString(a));
	}
	
	public static int getInt(){
		Scanner sc = new Scanner(System.in);
		int input = sc.nextInt();
		sc.close();
		return input;
	}
	
	public static double getDouble(){
		Scanner sc = new Scanner(System.in);
		double input = sc.nextDouble();
		sc.close();
		return input;
	}
	
	public static String getLine(){
		String input = sc.nextLine();
		return input;
	}
	
	public static String getNext(){
		String input = sc.next();
		return input;
	}
	
	public static String[] getList(){
		String input;
		List<String> sort = new ArrayList<String>(0);
		String StringBundle = "";
		char[] inputList;
		input = sc.nextLine();
		inputList = input.toCharArray();
		for(char a:inputList){
			if(a == ','){
				sort.add(StringBundle);
				StringBundle = "";
			}
			else StringBundle += Character.toString(a);
		}
		sort.add(StringBundle);
		return sort.toArray(new String[sort.size()]);
	}
	public static String[] getList(String Divider){
		String input;
		List<String> sort = new ArrayList<String>(0);
		String StringBundle = "";
		char[] inputList;
		input = sc.nextLine();
		inputList = input.toCharArray();
		for(char a:inputList){
			if(Divider.indexOf(a)>=0){
				if(StringBundle!="")
				sort.add(StringBundle);
				StringBundle = "";
			}
			else StringBundle += Character.toString(a);
		}
		if(StringBundle!="")
		sort.add(StringBundle);
		return sort.toArray(new String[sort.size()]);
	}
	public static String[] getList(char Divider){
		String input;
		List<String> sort = new ArrayList<String>(0);
		String StringBundle = "";
		char[] inputList;
		Scanner sc = new Scanner(System.in);
		input = sc.nextLine();
		inputList = input.toCharArray();
		for(char a:inputList){
			if(a == Divider){
				sort.add(StringBundle);
				StringBundle = "";
			}
			else StringBundle += Character.toString(a);
		}
		sort.add(StringBundle);
		sc.close();
		return sort.toArray(new String[sort.size()]);
	}
}