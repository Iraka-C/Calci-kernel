/**
 * Testing the package
 * Including main()
 */

public class TestMain{
	private final static String string="fzero(eval(1/2+1/(x-1)+2integ(sin(reg*arctan(x))/(1+x^2)^(reg/2)/(exp(2πx)-1),0,∞),reg(x)),14i)";
	public static void main(String[] args){
		Expression expr=new Expression(string);
		System.out.println(expr.value().val);
	}
}
