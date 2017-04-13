package complexExpression;

import android.graphics.Color;
import android.util.Log;

import com.iraka.calci.MainActivity;

import java.util.ArrayList;
import java.util.List;

public class Expression {
	public String text;
	public int[] br; // Bracket Depth
	public int[] lastLB; // last left bracket of the same level
	public int[] nextFS; // next functional symbol of the same level
	public int[] commaCnt; // for each '(', how many comma belongs to it?
	public int[] funcSer; // the function interpreted
	private int brDiff; // the difference of ( and )

	// Cache number parsing result, reduce time consumption
	private int[] numberParseRPos; // Where does a number starts at a position end? -1 for not a number start.
	private double[] numberParseResult; // The value of a number starts at a position

	// Cache symbol parsing result
	private int[] isAddSubParseResult; // is this an add/sub symbol ?
	private boolean[] isOmitMult; // is '*' omitted here ?

	public volatile boolean isWorking;
	public volatile boolean isExited;

	private static final String mathOperator="+-*/^=√";
	private static Complex memValue=new Complex(Double.NaN,Double.NaN); // for memory function

	/*
	In this instance, the following data structure is used:
	Find all functional symbols (L/R bracket, comma) quickly

	pos:  0123456789a
	text: (1,2,(3,4))
	next: 2 4 a7 9
	last: 0 0 05 5 50
	br:   01111122210
	ccnt: 2    1
	*/

	public Expression(String s){
		text=s;
		br=new int[s.length()+1];
		lastLB=new int[s.length()+1];
		nextFS=new int[s.length()+1];
		commaCnt=new int[s.length()+1];
		funcSer=new int[s.length()+1];
		brDiff=0;

		int[] symbolStack=new int[s.length()+1]; // a position stack of all left brackets
		int[] lastSymbol=new int[s.length()+1]; // what's the position of the last symbol ?

		int top=-1;

		br[0]=0;
		for(int i=0;i<s.length();i++){
			lastLB[i]=-1;
			nextFS[i]=-1;
			commaCnt[i]=0;
			funcSer[i]=-1;
			char c=s.charAt(i);
			if(i>0){
				br[i]=br[i-1];
				if(s.charAt(i-1)=='(')br[i]++;
				if(c==')')br[i]--;
			}

			if(c=='('){ // push
				top++;
				symbolStack[top]=i;
				lastLB[i]=i;
				lastSymbol[top]=i;
				brDiff++;
			}
			if(c==','&&top>=0){ // record
				lastLB[i]=symbolStack[top];
				commaCnt[symbolStack[top]]++;
				nextFS[lastSymbol[top]]=i;
				lastSymbol[top]=i;
			}
			if(c==')'&&top>=0){ // pop
				lastLB[i]=symbolStack[top];
				nextFS[lastSymbol[top]]=i;
				top--;
				brDiff--;
			}
		}

		numberParseRPos=new int[s.length()];
		numberParseResult=new double[s.length()];
		isAddSubParseResult=new int[s.length()];
		isOmitMult=new boolean[s.length()];

		for(int i=0;i<s.length();i++){
			numberParseRPos[i]=-1;
			isAddSubParseResult[i]=-1;
			isOmitMult[i]=false;
		}
	}

	private boolean isOperator(char c){
		return mathOperator.indexOf(c)>=0;
	}

	// only used when text[p]=='+'/'-' && p>0
	private boolean isAddSubSymbol(int p){
		if(isAddSubParseResult[p]>-1)
			return isAddSubParseResult[p]>0;

		char cj=text.charAt(p);
		if(!(cj=='+'||cj=='-')){
			isAddSubParseResult[p]=0;
			return false;
		}

		cj=text.charAt(p-1);
		if(isOperator(cj)||cj=='E'){
			isAddSubParseResult[p]=0;
			return false;
		}
		if(ParseNumber.isBaseSymbol(cj)){ // a pos/neg symbol in scientific notation under certain base
			int pos;
			for(pos=p+1;pos<text.length();pos++){
				cj=text.charAt(pos);
				if(!(cj>='0'&&cj<='9')){ // not a decimal number
					break;
				}
			}
			if(pos==text.length()){ // parsed to an end
				isAddSubParseResult[p]=0;
				return false;
			}
			if(pos==p+1){ // '+/-' directly followed by non-integer symbol
				isAddSubParseResult[p]=1;
				return true;
			}
			if(ParseNumber.isBaseSymbol(cj)||(cj>='A'&&cj<='F')||cj=='.'){ // part of another notation
				isAddSubParseResult[p]=1;
				return true;
			}
			isAddSubParseResult[p]=0;
			return false;
		}
		isAddSubParseResult[p]=1;
		return true;
	}

	// 0+NaN*I is never possible during a calculation
	// and is so used as "No Variable X provided" sign
	public Result value(int l,int r,Complex vX){
		if(!isWorking)return new Result(3); // Calculation Thread Halted

		if(l>r){
			return new Result(2).append(
			"Interpreter",
			"Expression expected",
			l,l
			);
		}

		String s=text.substring(l,r+1);
		// Constants
		if(s.equals("e"))return new Result(Complex.E); // constant e
		if(s.equals("pi")||s.equals("π"))return new Result(Complex.PI); // constant pi
		if(s.equals("i"))return new Result(Complex.I); // constant i
		if(s.equals("∞"))return new Result(Complex.Inf); // constant Infinity
		if(s.equals("inf"))return new Result(new Complex(Double.POSITIVE_INFINITY));
		if(s.equals("nan"))return new Result(new Complex(Double.NaN));
		if(s.equals("reg"))return new Result(memValue); // reg value
		if(s.equals("°"))return new Result(new Complex(Math.PI/180)); // degree value
		if(s.equals("%"))return new Result(new Complex(0.01)); // percent value
		if(s.equals("x")){
			if(vX.isValid()||vX.isNaN())
				return new Result(vX); // variable X
		}

		// try to parse s as a number
		if(numberParseRPos[l]==r){ // already parsed and cached
			return new Result(new Complex(numberParseResult[l]));
		}
		else try{
			// Forbid default real values and char e as a operator
			/*if(s.equals("Infinity")||s.equals("+Infinity")||s.equals("-Infinity")||
			s.equals("NaN")||s.equals("+NaN")||s.equals("-NaN")||s.indexOf('e')>=0){*/

			if(s.indexOf('e')>=0){
				throw new NumberFormatException();
			}

			try{ // try parse decimal double
				if(s.indexOf('D')>=0||s.indexOf('F')>=0){ // forbid double and float sign
					throw new NumberFormatException();
				}

				double v=Double.parseDouble(s);
				numberParseRPos[l]=r;
				numberParseResult[l]=v;
				return new Result(new Complex(v));
			}catch(NumberFormatException e){ // try parse double under a base
				// Not a valid dec Double
				double v=ParseNumber.parse(s);
				numberParseRPos[l]=r;
				numberParseResult[l]=v;
				return new Result(new Complex(v));
			}
		}catch(NumberFormatException e){
			// Not a valid Number
		}

		// omit space and enter
		if(text.charAt(l)==' '||text.charAt(l)=='\n'||text.charAt(l)=='\r'){
			return value(l+1,r,vX);
		}
		if(text.charAt(r)==' '||text.charAt(r)=='\n'||text.charAt(r)=='\r'){
			return value(l,r-1,vX);
		}

		char ci;
		char cj;

		// Addition and Subtraction
		for(int i=r;i>l;i--){
			ci=text.charAt(i);
			// Only ONE of the following long boolean expression will be calculated
			if(br[i]==br[l]&&isAddSubSymbol(i)){
				if(ci=='+'){
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.add(r1.val,r2.val));
				}
				else if(ci=='-'){
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.sub(r1.val,r2.val));
				}
			}
		}

		// Unary operator: positive and negative
		if(text.charAt(l)=='+')
			return value(l+1,r,vX);
		else if(text.charAt(l)=='-'){
			Result r1=value(l+1,r,vX);
			if(r1.isFatalError())return r1;
			return new Result(Complex.inv(r1.val));
		}

		// Multiplication and Division
		for(int i=r;i>l;i--){
			if(br[i]==br[l]){
				ci=text.charAt(i);
				cj=text.charAt(i-1);
				if(ci=='*'||ci=='×'){
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.mul(r1.val,r2.val));
				}
				else if(ci=='/'||ci=='÷'){
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.div(r1.val,r2.val));
				}
				else{ // * symbol omission
					if(isOmitMult[i]){ // cached
						Result r1=value(l,i-1,vX);
						if(r1.isFatalError()) return r1;
						Result r2=value(i,r,vX);
						if(r2.isFatalError()) return r2;
						return new Result(Complex.mul(r1.val,r2.val));
					}

					boolean iscjPreSymbol=(cj==')'||cj=='∞'||cj=='π'||cj=='°'||cj=='%');
					boolean iscjNumber=(cj>='0'&&cj<='9'||cj=='.');
					boolean iscjBase=ParseNumber.isBaseSymbol(cj);
					boolean case1=(ci>='a'&&ci<='z'||ci=='_'||ci=='(')&&(iscjNumber||iscjPreSymbol||iscjBase);
					boolean case2=(ci>='0'&&ci<='9'||ci=='.')&&(iscjPreSymbol);
					boolean case3=(ci=='∞'||ci=='π'||ci=='°'||ci=='%'||ci=='√'||ci=='Γ')&&(iscjNumber||iscjPreSymbol||iscjBase||cj>='a'&&cj<='z'||cj=='_');

					if(case1||case2||case3){
						isOmitMult[i]=true;
						Result r1=value(l,i-1,vX);
						if(r1.isFatalError()) return r1;
						Result r2=value(i,r,vX);
						if(r2.isFatalError()) return r2;
						return new Result(Complex.mul(r1.val,r2.val));
					}
				}
			}
		}

		// Power (priority right->left)
		for(int i=l;i<=r;i++)
			if(br[i]==br[l]&&text.charAt(i)=='^'){
				Result r1=value(l,i-1,vX);
				if(r1.isFatalError())return r1;
				Result r2=value(i+1,r,vX);
				if(r2.isFatalError())return r2;
				Complex t1=r1.val;
				Complex t2=r2.val;
				if(t1.re==0&&t1.im==0){ // special treatment for 0
					if(t2.re>0)return new Result(new Complex(0));
					else if(t2.re<0&&t2.im==0)return new Result(Complex.Inf);
					else return new Result(new Complex(Double.NaN,Double.NaN));
				}
				if(t1.norm().re<1&&t2.re==Double.POSITIVE_INFINITY){ // special treatment for inf
					return new Result(new Complex(0));
				}
				if(t1.norm().re>1&&t2.re==Double.NEGATIVE_INFINITY){ // special treatment for -inf
					return new Result(new Complex(0));
				}

				return new Result(Complex.exp(Complex.mul(t2,Complex.ln(t1))));
			}

		// Sqrt symbol
		if(text.charAt(l)=='√'){
			Result r1=value(l+1,r,vX);
			if(r1.isFatalError())return r1;
			return new Result(Complex.sqrt(r1.val));
		}

		// Brackets
		if(text.charAt(r)!=')')return new Result(1).append("Evaluator",s+" can't be calculated",l,r);
		if(text.charAt(l)=='(')return value(l+1,r-1,vX);

		// Functions
		int listPos; // the Position in Function class
		int funcID; // the ID in Function class
		int paramNum; // how many params in the function
		int leftBr; // where's the left bracket
		int exprParamNum; // how many params requires functional input

		if(funcSer[l]<0){ // not searched in list yet
			for(int i=0;i<Function.funcList.length;i++){
				if(s.startsWith(Function.funcList[i].funcName+"(")){
					funcSer[l]=i; // found
					break;
				}
			}
		}

		listPos=funcSer[l];

		// Not found
		if(listPos<0)return new Result(1).append("Evaluator",s+" can't be calculated",l,r);

		funcID=Function.funcList[listPos].funcSerial;
		leftBr=l+Function.funcList[listPos].funcName.length();
		exprParamNum=Function.funcList[listPos].exprParamNum;
		if(text.charAt(leftBr+1)==')'){
			paramNum=0;
		}
		else{
			paramNum=commaCnt[leftBr]+1;
		}

		//Log.i("Eval","s="+s+" l="+l+" len="+Function.funcList[i].funcName.length());
		//Log.i("Eval","ID="+funcID+" param="+paramNum+" leftbr="+leftBr+" exprp.="+exprParamNum);

		// Too many param.
		if(paramNum>9)return new Result(1).append("Interpreter","Too many parameters for function "+Function.funcList[listPos].funcName,l,r);

		// Calculate each param. value
		Complex[] val=new Complex[10];
		if(paramNum>0){
			for(int p=leftBr,i=0;nextFS[p]>=0;p=nextFS[p],i++){
				if(i>=exprParamNum){
					int resl=p+1;
					int resr=nextFS[p]-1;
					Result res=value(resl,resr,vX);
					if(res.isFatalError())
						return res.append("Evaluator","Invalid value for function "+Function.funcList[listPos].funcName,resl,resr);
					val[i]=res.val;
				}
			}
		}

		int funcJump=funcID+paramNum;
		switch(funcJump){
			case Function.EXP+1:return new Result(Complex.exp(val[0]));
			case Function.LN+1:return new Result(Complex.ln(val[0]));
			case Function.RE+1:return new Result(new Complex(val[0].re));
			case Function.IM+1:return new Result(new Complex(val[0].im));
			case Function.SQRT+1:return new Result(Complex.sqrt(val[0]));
			case Function.ABS+1:return new Result(val[0].norm()).append("Interpreter","Function abs is deprecated, please use norm instead",l,r);
			case Function.NORM+1:return new Result(val[0].norm());
			case Function.ARG+1:return new Result(val[0].arg());
			case Function.SIN+1:return new Result(Complex.sin(val[0]));
			case Function.COS+1:return new Result(Complex.cos(val[0]));
			case Function.TAN+1:return new Result(Complex.tan(val[0]));
			case Function.ARCSIN+1:return new Result(Complex.arcsin(val[0]));
			case Function.ARCCOS+1:return new Result(Complex.arccos(val[0]));
			case Function.ARCTAN+1:return new Result(Complex.arctan(val[0]));
			case Function.GAMMA+1:return new Result(Complex.gamma(val[0]));
			case Function.FLOOR+1:return new Result(new Complex(Math.floor(val[0].re),Math.floor(val[0].im)));
			case Function.CEIL+1:return new Result(new Complex(Math.ceil(val[0].re),Math.ceil(val[0].im)));
			case Function.REG:return new Result(memValue);
			case Function.REG+1:memValue=val[0];return new Result(val[0]);
			case Function.CONJ+1:return new Result(new Complex(val[0].re,-val[0].im));
			case Function.RAND:return new Result(new Complex(Math.random(),Math.random()));
			case Function.RAND+1:return new Result(new Complex(val[0].re*Math.random(),val[0].im*Math.random()));
			case Function.ROUND+1:return new Result(new Complex(Math.round(val[0].re),Math.round(val[0].im)));
			case Function.ROUND+2:
				double precRnd=Math.round(val[1].re);
				if(precRnd<0)return new Result(-1).append("Evaluator","Precision too low",l,r);
				if(precRnd>15)return new Result(-1).append("Evaluator","Precision too high",l,r);
				double ratio=Math.pow(10,precRnd);
				return new Result(new Complex(Math.round(val[0].re*ratio)/ratio,Math.round(val[0].im*ratio)/ratio));
			case Function.DIFF+2:return diff(leftBr+1,nextFS[leftBr]-1,val[1]);
			case Function.DIFF+3:return diff(leftBr+1,nextFS[leftBr]-1,val[1],val[2]);
			case Function.LIMIT+2:return limit(leftBr+1,nextFS[leftBr]-1,val[1]);
			case Function.LIMIT+3:return limit(leftBr+1,nextFS[leftBr]-1,val[1],val[2]);
			case Function.EVAL+2:return value(leftBr+1,nextFS[leftBr]-1,val[1]);
			case Function.FZERO+2:return solve(leftBr+1,nextFS[leftBr]-1,val[1]);
			case Function.INTEG+3:return integrate(leftBr+1,nextFS[leftBr]-1,val[1],val[2]);
			case Function.SUM+3:return sum(leftBr+1,nextFS[leftBr]-1,val[1],val[2]);
			case Function.PERM+1:return new Result(Complex.gamma(new Complex(val[0].re+1,val[0].im)));
			case Function.PERM+2:return new Result(perm(val[0],val[1]));
			case Function.COMB+2:return new Result(comb(val[0],val[1]));
			case Function.PREC:
				MainActivity.setGlobalPrecision(10);
				Result.setBase(Result.base);
				return new Result(new Complex(0)).append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			case Function.PREC+1:
				int prec=(int)Math.round(val[0].re);
				if(prec<0)return new Result(-1).setVal(new Complex(-1)).append("Setting","Precision too low",l,r);
				if(prec>Result.maxPrecision)
					return new Result(-1).setVal(new Complex(-1))
					.append("Setting","Precision too high",l,r)
					.append("Setting","The maximum precision is "+Result.maxPrecision+" digits",l,r);
				MainActivity.setGlobalPrecision(prec);
				Result.precision=prec;
				return new Result(new Complex(0)).append("Setting","Precision is set to "+prec+" digits",l,r);
			case Function.COLOR:
				int color=Color.parseColor("#74b1e2");
				MainActivity.setGlobalColor(color);
				return new Result(0).setVal(new Complex(0)).append("Interpreter","Reset color",l,r);
			case Function.COLOR+3:
				int vR=(int)val[0].re; // R
				if(vR<0||vR>255)return new Result(-1).append("Setting","Invalid red value provided",l,r);
				int vG=(int)val[1].re; // G
				if(vG<0||vG>255)return new Result(-1).append("Setting","Invalid green value provided",l,r);
				int vB=(int)val[2].re; // B
				if(vB<0||vB>255)return new Result(-1).append("Setting","Invalid blue value provided",l,r);

				int cMin=Math.min(vR,Math.min(vG,vB));
				int cMax=Math.max(vR,Math.max(vG,vB));
				double saturation=(cMax>0?1-(double)cMin/cMax:0);
				double brightness=cMax/255.0;

				if(saturation<0.3||brightness<0.5){
					return new Result(-1).setVal(new Complex(-1))
					.append("ColorPalette","Color assignment failed",l+6,r-1)
					.append("ColorPalette","This color may cause confusion",l,r);
				}

				color=Color.argb(255,vR,vG,vB);
				MainActivity.setGlobalColor(color);
				return new Result(0).setVal(new Complex(0)).append("Setting","Color is set to (R,G,B)=("+vR+","+vG+","+vB+")",l,r);
			case Function.BASE:
				Result.setBase(10);
				return new Result(new Complex(0))
				.append("Setting","Base is set to base-"+10+" ("+ParseNumber.baseName[10]+")",l,r)
				.append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			case Function.BASE+1:
				int base=(int)Math.round(val[0].re);
				if(!(base>=2&&base<=10||base==12||base==16))
					return new Result(-1).setVal(new Complex(-1)).append("Setting","Invalid base value",l,r);
				Result.setBase(base);
				return new Result(new Complex(0))
				.append("Setting","Base is set to base-"+base+" ("+ParseNumber.baseName[base]+")",l,r)
				.append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			default:
		}

		return new Result(1).append("Interpreter","Function "+Function.funcList[listPos].funcName+" syntax error",l,r);

	}

	public Result value(){ // Entrance !
		// 0+NaNi id a sign for "No Variable X provided" initially

		isWorking=true;
		isExited=false;

		if(brDiff!=0){
			isExited=true;
			return new Result(1).append("Interpreter","Brackets not paired",text.length()-1,text.length()-1);
		}
		Result res=value(0,text.length()-1,new Complex(0,Double.NaN));
		isExited=true;
		return res;
	}

	// 3 point diff()
	Result diff3(int l,int r,Complex x0,Complex delta){
		Result rn=value(l,r,Complex.sub(x0,delta));
		if(rn.isFatalError())return rn;
		if(!rn.val.isValid())return new Result(-1);
		Result rp=value(l,r,Complex.add(x0,delta));
		if(rp.isFatalError())return rp;

		Complex dv=Complex.div(Complex.sub(rp.val,rn.val),new Complex(delta.re*2,delta.im*2));
		return new Result(dv);
	}

	// 5 point diff()
	Result diff5(int l,int r,Complex x0,Complex delta){
		Result r1=diff3(l,r,x0,delta);
		if(r1.isFatalError())return r1;
		if(!r1.val.isValid())return new Result(-1);
		Result r2=diff3(l,r,x0,new Complex(delta.re*2,delta.im*2));
		if(r2.isFatalError())return r2;

		Complex dv=Complex.div(new Complex(r1.val.re*4-r2.val.re,r1.val.im*4-r2.val.im),new Complex(3));
		return new Result(dv);
	}

	// general diff()
	private boolean isDiffOverTolerance=false; // only checked for once, reduce data traffic
	Result diff(int l,int r,Complex x0){

		final int sect=8;
		final double sectAngle=Math.PI/sect;
		final double TOL=1E-5;
		Complex[] dirDer=new Complex[sect];
		Complex dsum=new Complex(0);
		double dvar=0;
		for(int i=0;i<sect;i++){ // find directions in [0,pi)
			Complex delta=new Complex(Math.cos(i*sectAngle)*TOL,Math.sin(i*sectAngle)*TOL);
			Result rv=diff5(l,r,x0,delta);
			if(rv.isFatalError())return rv;
			dirDer[i]=rv.val;
			dsum=Complex.add(dsum,rv.val);
		}
		dsum.re/=sect;
		dsum.im/=sect;
		for(int i=0;i<sect;i++){
			dvar+=Complex.sub(dirDer[i],dsum).norm2();
		}
		Result res=new Result(dsum);
		if(!isDiffOverTolerance&&dvar>TOL){
			isDiffOverTolerance=true;
			res.append("Differentiator","Function might be ill-conditioned at certain interval",l,r);
		}
		return res;
	}

	// directional diff
	Result diff(int l,int r,Complex x0,Complex dir){
		if(dir.re==0&&dir.im==0||!dir.isFinite())
			return new Result(1).append("Differentiator","Invalid direction value",l,r);

		final double TOL=1E-5;

		double norm=dir.norm().re;
		Complex delta=new Complex(dir.re/norm*TOL,dir.im/norm*TOL);
		Result rv=diff5(l,r,x0,delta);
		if(rv.isFatalError())return rv;
		return new Result(rv.val);
	}

	// solve under condition
	Result solve(int l,int r,Complex x0,Complex M,int iter){
		// Improved Secant Method, avoid second derivative
		Complex x1=x0;
		Result res1=value(l,r,x1);
		if(res1.isFatalError())return res1;
		Complex v1=res1.val;
		Complex r1=Complex.div(v1,diff(l,r,x1).val);
		if(r1.isNaN()){
			return new Result(-1).append("Solver","Invalid initial value",l,r); // Error occurred
		}
		// 1 step Newton Method Iteration
		Complex x2=Complex.sub(x1,r1);
		Complex v2=value(l,r,x2).val;
		Complex r2=Complex.div(v2,diff(l,r,x2).val);
		if(r2.isNaN()){
			return new Result(-1).append("Solver","Invalid initial value",l,r); // Error occurred
		}
		if(Complex.sub(x1,x2).norm2()<1E-20&&v2.norm2()<1E-20){ // 1E-10 precision
			return new Result(x2);
		}
		Complex x3;

		Result root=new Result(0);

		for(int i=0;i<=iter;i++){ // normally no more than 20 iter., but for eq. such as x^.2, more is needed.

			if(!isWorking)return new Result(3); // Calculation Thread Halted

			Complex d1=Complex.mul(Complex.sub(x2,x1),r2);
			Complex d2=Complex.sub(r2,r1);

			x3=Complex.mul(Complex.sub(x2,Complex.div(d1,d2)),M); // Relaxation Param.
			v2=value(l,r,x3).val;

			if(Complex.sub(x2,x3).norm2()<1E-20){ // 1E-10 precision
				if(v2.norm2()<1E-20){
					return new Result(x3);
				}
			}

			x1=x2;x2=x3;
			r1=r2;r2=Complex.div(v2,diff(l,r,x3).val);
			if(r2.isNaN()){ // Math Error
				if(M.re==1.0){ // 1.0 is a safely expressed double
					root=new Result(-1)
					.append("Solver","Invalid iteration occurred",l,r)
					.append("Solver","Function might be ill-conditioned at certain interval",l,r);
				}
				break;
			}
		}

		return root; // Not Found
	}

	// general solve
	Result solve(int l,int r,Complex x0){ // auto solver
		Result rp;
		for(double M=1.0;M>0.05;M*=0.7){

			if(!isWorking)return new Result(3); // Calculation Thread Halted
			//System.out.println("Relaxation M = "+new Double(M).toString()
			//		+" , MaxLoop = "+new Integer((int)Math.round(1500/Math.sqrt(M))).toString());
			rp=solve(l,r,x0,new Complex(M),(int)Math.round(1500/Math.sqrt(M)));
			if(rp.isFatalError())return rp;
			if(rp.val.isValid()&&!rp.val.isNaN()){
				return rp;
			}
			if(rp.err==-1){ // Initial Value Error
				break;
			}
			if(M==1.0){
				rp.append("Solver","Try Under-Relaxation Method",l,r);
			}
		}
		rp=new Result(-1);
		rp.append("Solver","Finding zero point "+text.substring(l,r+1)+" failed",l,r);
		return rp;
	}

	// simpson's integration from x0 to x2
	Complex simpsonIntegrate(int l,int r,Complex x0,Complex x2){
		Complex length=Complex.sub(x2,x0);
		Complex x1=new Complex((x0.re+x2.re)/2,(x0.im+x2.im)/2);
		Complex y0=value(l,r,x0).val;
		Complex y1=value(l,r,x1).val;
		Complex y2=value(l,r,x2).val;
		Complex sum=new Complex((y0.re+4*y1.re+y2.re)/6,(y0.im+4*y1.im+y2.im)/6);
		sum=Complex.mul(sum,length);

		return sum;
	}

	// gauss's integration from x0 to x2
	private static final double giIntv=Math.sqrt(0.6);
	Complex gaussIntegrate3(int l,int r,Complex x0,Complex x2){
		Complex length=Complex.sub(x2,x0);
		Complex x1n=new Complex((x0.re+x2.re)/2,(x0.im+x2.im)/2);
		Complex x0n=new Complex(x1n.re-giIntv*length.re/2,x1n.im-giIntv*length.im/2);
		Complex x2n=new Complex(x1n.re+giIntv*length.re/2,x1n.im+giIntv*length.im/2);
		Complex y0=value(l,r,x0n).val;
		Complex y1=value(l,r,x1n).val;
		Complex y2=value(l,r,x2n).val;
		Complex sum=new Complex((5*y0.re+8*y1.re+5*y2.re)/18,(5*y0.im+8*y1.im+5*y2.im)/18);
		sum=Complex.mul(sum,length);

		return sum;
	}

	// gauss integration at 5 points
	private static final double gi5Intv1=Math.sqrt(5-2*Math.sqrt(10./7.))/3;
	private static final double gi5Intv2=Math.sqrt(5+2*Math.sqrt(10./7.))/3;
	private static final double gi5K0=128./225.;
	private static final double gi5K1=(322+13*Math.sqrt(70))/900;
	private static final double gi5K2=(322-13*Math.sqrt(70))/900;
	Complex gaussIntegrate(int l,int r,Complex x0,Complex x2){
		Complex lenH=new Complex(x2.re-x0.re,x2.im-x0.im);
		Complex halfH=new Complex(lenH.re/2,lenH.im/2);

		Complex t1=new Complex((x0.re+x2.re)/2,(x0.im+x2.im)/2);
		Complex t2=new Complex(t1.re+gi5Intv1*halfH.re,t1.im+gi5Intv1*halfH.im);
		Complex t3=new Complex(t1.re+gi5Intv2*halfH.re,t1.im+gi5Intv2*halfH.im);
		Complex t4=new Complex(t1.re-gi5Intv1*halfH.re,t1.im-gi5Intv1*halfH.im);
		Complex t5=new Complex(t1.re-gi5Intv2*halfH.re,t1.im-gi5Intv2*halfH.im);

		Complex y1=value(l,r,t1).val;
		Complex y2=value(l,r,t2).val;
		Complex y3=value(l,r,t3).val;
		Complex y4=value(l,r,t4).val;
		Complex y5=value(l,r,t5).val;

		Complex sum=new Complex(
		y1.re*gi5K0+(y2.re+y4.re)*gi5K1+(y3.re+y5.re)*gi5K2,
		y1.im*gi5K0+(y2.im+y4.im)*gi5K1+(y3.im+y5.im)*gi5K2
		);
		sum=Complex.mul(sum,halfH);

		return sum;
	}

	// integrate from x0 to x2, auto step-length
	private boolean isIntegOverTolerance=false; // only checked for once, reduce data traffic
	Result adaptiveIntegrate(int l,int r,Complex x0,Complex x2,Complex lastSum,double TOL,int depth){

		if(!isWorking)return new Result(3); // Calculation Thread Halted

		Complex x1;

		if(Double.isInfinite(x0.re)){
			return new Result(3).setVal(new Complex(0));
		}

		if(Double.isInfinite(x2.re)){
			double aRe=Math.abs(x0.re);
			double newRe;
			if(aRe<=5000) // Avoid precision lost due to large number
				newRe=(aRe<1?Math.exp(aRe):Math.E*aRe);
			else
				return new Result(-1).setVal(new Complex(0));

			//System.out.println("New Point: "+x0.toString()+" ~ "+x2.toString());
			if(x2.re>0){ // POS Inf
				x1=new Complex(newRe,x0.im);
			}
			else{ // NEG Inf
				x1=new Complex(-newRe,x0.im);
			}
		}
		else{
			x1=new Complex((x0.re+x2.re)/2,(x0.im+x2.im)/2);
		}

		Complex sAB,sAC,sCB,sABnew,abbr;

		sAB=lastSum;
		sAC=gaussIntegrate(l,r,x0,x1);
		sCB=gaussIntegrate(l,r,x1,x2);
		sABnew=new Complex(0);
		if(sAC.isFinite())sABnew=Complex.add(sABnew,sAC);
		if(sCB.isFinite())sABnew=Complex.add(sABnew,sCB);
		abbr=Complex.sub(sAB,sABnew);

		if(abbr.isValid()&&abbr.norm2()<200*TOL){
			return new Result(sABnew);
		}

		if(depth>=18){
			Result r1=new Result(sABnew);
			if(!isIntegOverTolerance&&abbr.norm2()>2E6*TOL){
				isIntegOverTolerance=true;
				r1.append("Quadrature","Function might be ill-conditioned at certain interval",l,r);
			}

			return r1;
		}

		sAC=adaptiveIntegrate(l,r,x0,x1,sAC,TOL/4,depth+1).val;
		sCB=adaptiveIntegrate(l,r,x1,x2,sCB,TOL/4,depth+1).val;
		sABnew=Complex.add(sAC,sCB);
		return new Result(sABnew);
	}

	// general quadrature
	Result integrate(int l,int r,Complex x0,Complex x2){
		Result check=value(l,r,x0);
		if(check.isFatalError())return check;

		if(Double.isInfinite(x0.re)){
			if(Double.isInfinite(x2.re)){ // Add middle point
				double iim=((Complex.isDoubleFinite(x0.im)?x0.im:0)+(Complex.isDoubleFinite(x2.im)?x2.im:0))/2;

				Result r1=integrate(l,r,x0,new Complex(0,iim));
				Result r2=integrate(l,r,new Complex(0,iim),x2);
				return new Result(Complex.add(r1.val,r2.val));
			}else{
				Result r1=integrate(l,r,x2,x0);
				return new Result(new Complex(-r1.val.re,-r1.val.im));
			}
		}

		double TOL=1E-8; // the precision expected, 8 digits

		Complex sAB=gaussIntegrate(l,r,x0,x2);
		return adaptiveIntegrate(l,r,x0,x2,sAB,TOL*TOL,0);
	}

	// summation
	Result sum(int l,int r,Complex start,Complex end){
		double ds=start.re;
		double de=end.re;
		if(de<ds){
			return new Result(2).append("Summation","Upper bound is smaller than lower bound",l,r);
		}

		Complex sum=new Complex(0);
		Complex v=new Complex(0);
		final double TOL2=1E-14; // 7 digits expected
		final int maxBoundCnt=1000;
		final int maxCnt=100000;
		int boundCnt=0;
		int cnt=0;

		double ratio=(end.im-start.im)/(de-ds);
		if(!Complex.isDoubleFinite(ratio)){
			return new Result(3).append("Summation","Cannot decide sum path",l,r);
		}
		for(v.re=ds;v.re<=de;v.re+=1,cnt++){
			if(!isWorking)return new Result(3); // Calculation Thread Halted

			v.im=(v.re-ds)*ratio+start.im;
			Result res=value(l,r,v);
			if(res.isFatalError()){
				return res;
			}
			if(!res.val.isFinite()){
				return new Result(sum)
				.append("Summation","Error occurred during summation at x="+v.toString(),l,r)
				.append("Summation","The sum might not be finite",l,r);
			}
			if(res.val.norm2()<TOL2){
				boundCnt++;
			}
			else{
				boundCnt=0;
			}
			if(boundCnt>maxBoundCnt){
				break;
			}
			if(cnt==maxCnt){
				new Result(-1).append("Summation","Time consumption is longer than expected",l,r);
			}

			sum=Complex.add(sum,res.val);
		}

		return new Result(sum);
	}

	// permutation & combination
	private Complex permIter(Complex n_,Complex m_){ // Gamma(n+1)/Gamma(m+1)
		Complex n,m;
		Complex ans=new Complex(1);
		n=n_;
		m=m_;

		for(;;){ // adapted from iteration
			if(n.re>1&&m.re>1){
				if(n.re-m.re>=1){
					ans=Complex.mul(new Complex(n.re),ans);
					n.re-=1;
				}
				else if(m.re-n.re>=1){
					ans=Complex.div(ans,new Complex(m.re));
					m.re-=1;
				}
				else{
					ans=Complex.mul(new Complex(n.re/m.re),ans);
					n.re-=1;
					m.re-=1;
				}
			}
			else if(n.re==m.re&&n.im==m.im){
				break;
			}
			else{
				Complex af=Complex.div(Complex.gamma(new Complex(n.re+1,n.im)),Complex.gamma(new Complex(m.re+1,m.im)));
				ans=Complex.mul(af,ans);
				break;
			}
			if(!ans.isFinite()){ // invalid value occurred, no need to continue
				break;
			}
			if(!isWorking){
				return new Complex(); // Calculation Thread Halted
			}
		}

		return ans;
	}
	public Complex perm(Complex n,Complex r){

		return permIter(n,Complex.sub(n,r));
	}
	private Complex combIter(Complex n_,Complex m_){ // Gamma(n+1)/Gamma(m+1)/Gamma(n-m+1)

		Complex n,m;
		Complex ans=new Complex(1);
		n=n_;
		m=m_;

		for(;;){ // adapted from iteration
			if(n.re>1&&m.re>1){
				if(n.re-m.re>=1){
					ans=Complex.mul(new Complex(n.re),ans);
					n.re-=1;
				}
				else if(m.re-n.re>=1){
					ans=Complex.div(ans,new Complex(m.re));
					m.re-=1;
				}
				else{
					ans=Complex.mul(new Complex(n.re/m.re),ans);
					n.re-=1;
					m.re-=1;
				}
			}
			else{
				Complex af=Complex.div(perm(n,m),Complex.gamma(new Complex(m.re+1,m.im)));
				ans=Complex.mul(af,ans);
				break;
			}
			if(!ans.isFinite()){ // invalid value occurred, no need to continue
				break;
			}
			if(!isWorking){
				return new Complex(); // Calculation Thread Halted
			}
		}

		return ans;
	}
	public Complex comb(Complex n,Complex r){
		return combIter(n,r);
	}

	// limit (Lagrange interpolation)
	private static final Complex par2p=new Complex((1+Math.sqrt(2))/4);
	private static final Complex par2n=new Complex((1-Math.sqrt(2))/4);
	private static final Complex hRatio=new Complex(1+Math.sqrt(2));
	private Result limitH(int l,int r,Complex x0,Complex h){

		Complex f0;

		if(x0.isFinite()){

			Complex x1=Complex.add(x0,h);
			Complex x2=Complex.add(x0,Complex.mul(h,hRatio));
			Complex x3=Complex.sub(x0,h);
			Complex x4=Complex.sub(x0,Complex.mul(h,hRatio));

			Result r1=value(l,r,x1);
			Result r2=value(l,r,x2);
			Result r3=value(l,r,x3);
			Result r4=value(l,r,x4);

			if(r1.isFatalError())return r1;
			if(r2.isFatalError())return r2;
			if(r3.isFatalError())return r3;
			if(r4.isFatalError())return r4;

			Complex f1=r1.val;
			Complex f2=r2.val;
			Complex f3=r3.val;
			Complex f4=r4.val;

			Complex f13=Complex.mul(Complex.add(f1,f3),par2p);
			Complex f24=Complex.mul(Complex.add(f2,f4),par2n);
			f0=Complex.add(f13,f24);
		}
		else{
			double norm2=h.norm2();
			Complex N=new Complex(h.re/norm2,h.im/norm2);
			return value(l,r,N);
		}

		return new Result(f0);
	}

	// general limit
	public Result limit(int l,int r,Complex x0){

		//new Result(-1).append("Evaluator","Function limit() still requires improvement. Results are only for development.",l,r);

		// Responsibility recharge
		if(Double.isInfinite(x0.re)){ // to real infinity
			return limit(l,r,x0,new Complex(x0.re>0?1:-1));
		}
		if(Double.isInfinite(x0.im)){ // to imaginary infinity
			return limit(l,r,x0,new Complex(0,x0.im>0?1:-1));
		}

		final int sect=8;
		final double sectAngle=Math.PI/sect;
		Complex[] limitRes=new Complex[sect];
		Complex limitSum=new Complex(0);
		double limitVar=0;
		int validSect=0;

		for(int i=0;i<sect;i++){ // find directions in [0,pi)

			List<Complex> histRes=new ArrayList<>();
			double minDe=1E200;
			int minPos=-1;

			int cnt=0;
			int overErrorRangeCount=0;
			double h;
			for(h=1E-1;h>=1E-10;h*=0.9,cnt++){
				Complex delta=new Complex(Math.cos(i*sectAngle)*h,Math.sin(i*sectAngle)*h);
				Result resR=limitH(l,r,x0,delta);
				if(resR.isFatalError())return resR.append("Evaluator","Limit not found",l,r);
				Complex res=resR.val;

				if(cnt>0){
					double e=Complex.sub(res,histRes.get(cnt-1)).norm().re;

					if(e<minDe){
						overErrorRangeCount=0;
						minDe=e;
						minPos=cnt;
					}
					else{ // e>=minDe or e is NaN!
						overErrorRangeCount++;
					}
					if(overErrorRangeCount>20){
						break;
					}

					//System.out.println("h="+Double.toString(h)+" err="+e);
				}

				histRes.add(res);
			}

			if(minDe>1E-5){ // didn't found ?
				new Result(-1).append("Evaluator","Function may be not convergent in all directions",l,r);
			}else{ // found
				Complex minRes=histRes.get(minPos-1);
				limitSum=Complex.add(limitSum,minRes);
				limitRes[validSect]=minRes;
				validSect++;
			}
		}

		if(validSect==0)return new Result(-1).append("Evaluator","Limit not found",l,r);

		limitSum.re/=validSect;
		limitSum.im/=validSect;
		for(int i=0;i<validSect;i++){
			limitVar+=Complex.sub(limitRes[i],limitSum).norm2();
		}

		//Log.i("Limit","Dvar="+limitVar);

		Result res=new Result(limitSum);
		if(limitVar>1E-5){
			res.append("Evaluator","Function may have not converged to expected precision",l,r);
		}

		return res;
	}
	// directional limit
	public Result limit(int l,int r,Complex x0,Complex dir){
		if(dir.re==0&&dir.im==0||!dir.isFinite())
			return new Result(1).append("Differentiator","Invalid direction value",l,r);


		List<Complex> histRes=new ArrayList<>();
		double minDe=1E200;
		int minPos=-1;
		double norm=dir.norm().re;

		int cnt=0;
		int overErrorRangeCount=0;
		double h;
		for(h=1E-1;h>=1E-10;h*=0.9,cnt++){
			Complex delta=new Complex(dir.re/norm*h,dir.im/norm*h);
			Result resR=limitH(l,r,x0,delta);
			if(resR.isFatalError())return resR.append("Evaluator","Limit not found",l,r);
			Complex res=resR.val;

			if(cnt>0){
				double e=Complex.sub(res,histRes.get(cnt-1)).norm().re;

				if(e<minDe){
					overErrorRangeCount=0;
					minDe=e;
					minPos=cnt;
				}
				else{ // e>=minDe or e is NaN!
					overErrorRangeCount++;
				}
				if(overErrorRangeCount>20){
					break;
				}

				//System.out.println("h="+Double.toString(h)+" val="+res+" err="+e);
			}

			histRes.add(res);
		}

		Complex minRes=histRes.get(minPos-1);
		//Log.i("Limit","err="+minDe);

		if(minDe>1E-5){ // didn't found ?
			return new Result(minRes).append("Evaluator","Function may be not convergent at given point",l,r);
		}else{ // found
			return new Result(minRes);
		}
	}

	// stop evaluation process, can be called asynchronizely
	public void stopEvaluation(){
		isWorking=false;
	}
}
