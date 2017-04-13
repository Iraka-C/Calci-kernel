package complexExpression;

import java.util.*;

// Containing a complex result
public class Result{
	public Complex val;
	public int err;
	public static int precision=10;
	public static int base=10;
	public static int maxPrecision=15;
	/* Error No.:
	 * 0: No Error
	 * 1: Grammatical Error
	 * 2: Item Expected Error
	 * 3: Runtime Error
	 * -1: Mathematical Error
	 * <=-2: Special Functions
	 */

	public Result(Complex v){
		val=v;
		err=0;
	}
	public Result(int err_){
		if(err_<-1){ // special functions
			// do sth for special functions
			return;
		}
		val=new Complex(Double.NaN,Double.NaN);
		err=err_;
	}
	public Result append(String name,String text,int l,int r){
		// Well, best to be static, but I don't want to fix it
		System.out.println(name+":: "+text);
		return this;
	}
	public Result setVal(Complex v_){
		val=v_;
		return this;
	}

	public static void setBase(int base_){
		base=base_;
		precision=(int)Math.floor(35*Math.log(2)/Math.log(base_));
		maxPrecision=(int)Math.floor(52*Math.log(2)/Math.log(base_));
	}

	boolean isFatalError(){
		return err>0;
	}
}
