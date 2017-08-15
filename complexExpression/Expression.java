import java.util.ArrayList;
import java.util.List;

public class Expression {
	public String text;
	public int[] br; // Bracket Depth
	public int[] lastLB; // last left bracket of the same level
	public int[] nextFS; // next functional symbol of the same level
	public int[] commaCnt; // for each '(', how many comma belongs to it?
	public int[] funcSer; // the function interpreted
	public int brDiff; // the difference of ( and )

	// Cache interpretation result
	private class SymbolCachePair{
		static final int SYMBOL_NUM=0; // may cache value
		static final int SYMBOL_ADD=1;
		static final int SYMBOL_POS=2;
		static final int SYMBOL_SUB=3;
		static final int SYMBOL_NEG=4;
		static final int SYMBOL_MUL=5;
		static final int SYMBOL_DIV=6;
		static final int SYMBOL_MUL_OMIT=7;
		static final int SYMBOL_POW=8;
		static final int SYMBOL_SQRT=9;
		static final int SYMBOL_CONST=10; // may cache value
		static final int SYMBOL_FUNC=11;
		static final int SYMBOL_VAR=12; // Temporarily not used
		static final int SYMBOL_BRACKET=13;

		int end_pos;
		int symbol;
		int symbol_pos;
		Complex cachedValue;
		SymbolCachePair(int end_pos_,int symbol_,int symbol_pos_,Complex cachedValue_){
			end_pos=end_pos_;
			symbol=symbol_;
			symbol_pos=symbol_pos_;
			cachedValue=cachedValue_;
		}
	}
	private class SymbolCache{
		List<SymbolCachePair> list;
		SymbolCache(){
			list=new ArrayList<>();
		}
		void submit(int end_pos_,int symbol_,int symbol_pos_){
			list.add(new SymbolCachePair(end_pos_,symbol_,symbol_pos_,new Complex()));
			//Log.i("expression","Submit pos="+symbol_pos_+" type="+symbol_);
		}
		void submit(int end_pos_,int symbol_,Complex cachedValue_){
			list.add(new SymbolCachePair(end_pos_,symbol_,-1,cachedValue_));
			//Log.i("expression","Submit val="+cachedValue_+" type="+symbol_);
		}
		SymbolCachePair checkCache(int end_pos){
			for(int i=0;i<list.size();i++){
				SymbolCachePair pair=list.get(i);
				if(pair.end_pos==end_pos){
					return pair;
				}
			}
			return null;
		}
	}
	private SymbolCache[] interpretResult;

	public volatile boolean isWorking;
	public volatile boolean isExited;

	private static final String mathOperator="+-*/^=√";
	private static Complex memValue=new Complex(); // for memory function
	private static Complex ansValue=new Complex(); // for last answer

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

	public Expression(String s){ // primitive analysis to the expression
		text=s;
		br=new int[s.length()+1];
		lastLB=new int[s.length()+1];
		nextFS=new int[s.length()+1];
		commaCnt=new int[s.length()+1];
		brDiff=0;

		int[] symbolStack=new int[s.length()+1]; // a position stack of all left brackets
		int[] lastSymbol=new int[s.length()+1]; // what's the position of the last symbol ?

		int top=-1;

		br[0]=0;
		for(int i=0;i<s.length();i++){
			lastLB[i]=-1;
			nextFS[i]=-1;
			commaCnt[i]=0;
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
			if(c==')'){
				if(top>=0){ // pop
					lastLB[i]=symbolStack[top];
					nextFS[lastSymbol[top]]=i;
					top--;
				}
				brDiff--;
			}
		}
	}

	private void initCache(){ // init cache for evaluation
		funcSer=new int[text.length()];
		interpretResult=new SymbolCache[text.length()];
		for(int i=0;i<interpretResult.length;i++){
			interpretResult[i]=new SymbolCache();
			funcSer[i]=-1;
		}
	}

	private boolean isOperator(char c){
		return mathOperator.indexOf(c)>=0;
	}

	// only used when text[p]=='+'/'-' && p>0
	private boolean isAddSubSymbol(int p){
		if(p==0)return false;

		char cj=text.charAt(p);
		if(!(cj=='+'||cj=='-')){
			return false;
		}

		cj=text.charAt(p-1);
		if(isOperator(cj)||cj=='E'){
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
				return false;
			}
			if(pos==p+1){ // '+/-' directly followed by non-integer symbol
				return true;
			}
			if(ParseNumber.isBaseSymbol(cj)||(cj>='A'&&cj<='F')||cj=='.'){ // part of another notation
				return true;
			}
			return false;
		}
		return true;
	}

	// only used when text[p]=='*' && p>0
	private boolean isOmitMult(int p){
		if(p==0)return false;

		char ci=text.charAt(p);
		char cj=text.charAt(p-1);

		boolean iscjPreSymbol=(cj==')'||cj=='∞'||cj=='π'||cj=='°'||cj=='%');
		boolean iscjNumber=(cj>='0'&&cj<='9'||cj>='A'&&cj<='Z'||cj=='.');
		boolean iscjBase=ParseNumber.isBaseSymbol(cj);
		boolean iscjFunc=(cj>='a'&&cj<='z'||cj=='_');
		boolean isciNumber=(ci>='0'&&ci<='9'||ci>='A'&&ci<='Z'||ci=='.');
		//boolean isciBase=ParseNumber.isBaseSymbol(ci);

		boolean case1=(ci>='a'&&ci<='z'||ci=='_'||ci=='(')&&(iscjNumber||iscjPreSymbol||iscjBase);
		boolean case2=(isciNumber)&&(iscjPreSymbol||iscjFunc);
		boolean case3=(ci=='∞'||ci=='π'||ci=='°'||ci=='%'||ci=='√'||ci=='Γ')&&(iscjNumber||iscjPreSymbol||iscjBase||iscjFunc);

		return case1||case2||case3;
	}

	// 0+NaN*I is never possible during a calculation
	// and is so used as "No Variable X provided" sign
	private Result value(int l,int r,Complex vX){
		if(!isWorking)return new Result(3); // Calculation Thread Halted

		if(l>r){ // item expected
			return new Result(2).append(
			"Interpreter",
			"Expression expected",
			l,l
			);
		}

		// Check if result cached
		SymbolCachePair pair=interpretResult[l].checkCache(r);
		if(pair!=null){ // cached result
			Result r1,r2;
			switch(pair.symbol){ // No fatal error now
				case SymbolCachePair.SYMBOL_CONST:
				case SymbolCachePair.SYMBOL_NUM:
					return new Result(pair.cachedValue);
				case SymbolCachePair.SYMBOL_ADD:
					r1=value(l,pair.symbol_pos-1,vX);
					r2=value(pair.symbol_pos+1,r,vX);
					return new Result(Complex.add(r1.val,r2.val));
				case SymbolCachePair.SYMBOL_SUB:
					r1=value(l,pair.symbol_pos-1,vX);
					r2=value(pair.symbol_pos+1,r,vX);
					return new Result(Complex.sub(r1.val,r2.val));
				case SymbolCachePair.SYMBOL_POS:
					return value(l+1,r,vX);
				case SymbolCachePair.SYMBOL_NEG:
					r1=value(l+1,r,vX);
					return new Result(Complex.inv(r1.val));
				case SymbolCachePair.SYMBOL_MUL:
					r1=value(l,pair.symbol_pos-1,vX);
					r2=value(pair.symbol_pos+1,r,vX);
					return new Result(Complex.mul(r1.val,r2.val));
				case SymbolCachePair.SYMBOL_DIV:
					r1=value(l,pair.symbol_pos-1,vX);
					r2=value(pair.symbol_pos+1,r,vX);
					return new Result(Complex.div(r1.val,r2.val));
				case SymbolCachePair.SYMBOL_MUL_OMIT:
					r1=value(l,pair.symbol_pos-1,vX);
					r2=value(pair.symbol_pos,r,vX); // Attention to the pos!
					return new Result(Complex.mul(r1.val,r2.val));
				case SymbolCachePair.SYMBOL_POW:
					r1=value(l,pair.symbol_pos-1,vX);
					r2=value(pair.symbol_pos+1,r,vX);
					return new Result(Complex.pow(r1.val,r2.val));
				case SymbolCachePair.SYMBOL_SQRT:
					r1=value(l+1,r,vX);
					return new Result(Complex.sqrt(r1.val));
				case SymbolCachePair.SYMBOL_FUNC:
					return funcValue(l,r,vX);
				case SymbolCachePair.SYMBOL_BRACKET:
					return value(l+1,r-1,vX);
			}
		}

		// Interpret expression
		String s=text.substring(l,r+1);

		// Variable
		if(s.equals("x")&&(vX.isValid()||vX.isNaN()))return new Result(vX); // variable X
		if(s.equals("reg"))return new Result(memValue); // reg value

		// omit space and enter
		if(text.charAt(l)==' '||text.charAt(l)=='\n'||text.charAt(l)=='\r')return value(l+1,r,vX);
		if(text.charAt(r)==' '||text.charAt(r)=='\n'||text.charAt(r)=='\r')return value(l,r-1,vX);

		//======================= Below this line, string will be parsed only once ========================

		{ // Constants
			Complex complexConst=null;
			if(s.equals("e")) complexConst=Complex.E; // constant e
			if(s.equals("pi")||s.equals("π")) complexConst=Complex.PI; // constant pi
			if(s.equals("i")) complexConst=Complex.I; // constant i
			if(s.equals("∞")) complexConst=Complex.Inf; // constant Infinity
			if(s.equals("inf")) complexConst=new Complex(Double.POSITIVE_INFINITY);
			if(s.equals("nan")) complexConst=new Complex();
			if(s.equals("°")) complexConst=new Complex(Math.PI/180); // degree value
			if(s.equals("%")) complexConst=new Complex(0.01); // percent value
			if(s.equals("ans")) complexConst=ansValue; // ans is stable during 1 calculation
			if(complexConst!=null){
				interpretResult[l].submit(r,SymbolCachePair.SYMBOL_CONST,complexConst);
				return new Result(complexConst);
			}
		}

		// Number parsing
		try{
			// Forbid default real values and char e as a operator
			// Forbid default notations
			if(s.indexOf('e')>=0||s.indexOf('I')>=0||s.indexOf('N')>=0||
				s.indexOf('X')>=0||s.indexOf('P')>=0||s.indexOf('x')>=0||s.indexOf('p')>=0){
				throw new NumberFormatException();
			}

			try{ // try parse decimal double
				if(s.indexOf('D')>=0||s.indexOf('F')>=0){ // forbid double and float sign
					throw new NumberFormatException();
				}

				double v=Double.parseDouble(s);
				interpretResult[l].submit(r,SymbolCachePair.SYMBOL_NUM,new Complex(v));
				return new Result(new Complex(v));
			}catch(NumberFormatException e){ // try parse double under a base
				// Not a valid dec Double
				double v=ParseNumber.parse(s);
				interpretResult[l].submit(r,SymbolCachePair.SYMBOL_NUM,new Complex(v));
				return new Result(new Complex(v));
			}
		}catch(NumberFormatException e){
			// Not a valid Number
		}

		char ci;
		// Addition and Subtraction
		for(int i=r;i>l;i--){
			ci=text.charAt(i);
			// Only ONE of the following long boolean expression will be calculated
			if(br[i]==br[l]&&isAddSubSymbol(i)){
				if(ci=='+'){
					interpretResult[l].submit(r,SymbolCachePair.SYMBOL_ADD,i);
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.add(r1.val,r2.val));
				}
				else if(ci=='-'){
					interpretResult[l].submit(r,SymbolCachePair.SYMBOL_SUB,i);
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.sub(r1.val,r2.val));
				}
			}
		}

		// Unary operator: positive and negative
		if(text.charAt(l)=='+'){
			interpretResult[l].submit(r,SymbolCachePair.SYMBOL_POS,-1);
			return value(l+1,r,vX);
		}
		else if(text.charAt(l)=='-'){
			interpretResult[l].submit(r,SymbolCachePair.SYMBOL_NEG,-1);
			Result r1=value(l+1,r,vX);
			if(r1.isFatalError())return r1;
			return new Result(Complex.inv(r1.val));
		}

		// Multiplication and Division
		for(int i=r;i>l;i--){
			if(br[i]==br[l]){
				ci=text.charAt(i);
				if(ci=='*'||ci=='×'){
					interpretResult[l].submit(r,SymbolCachePair.SYMBOL_MUL,i);
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.mul(r1.val,r2.val));
				}
				else if(ci=='/'||ci=='÷'){
					interpretResult[l].submit(r,SymbolCachePair.SYMBOL_DIV,i);
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError())return r1;
					Result r2=value(i+1,r,vX);
					if(r2.isFatalError())return r2;
					return new Result(Complex.div(r1.val,r2.val));
				}
				else if(isOmitMult(i)){ // * symbol omission
					interpretResult[l].submit(r,SymbolCachePair.SYMBOL_MUL_OMIT,i);
					Result r1=value(l,i-1,vX);
					if(r1.isFatalError()) return r1;
					Result r2=value(i,r,vX);
					if(r2.isFatalError()) return r2;
					return new Result(Complex.mul(r1.val,r2.val));
				}
			}
		}

		// Power (priority right->left)
		for(int i=l;i<=r;i++)
			if(br[i]==br[l]&&text.charAt(i)=='^'){
				interpretResult[l].submit(r,SymbolCachePair.SYMBOL_POW,i);
				Result r1=value(l,i-1,vX);
				if(r1.isFatalError())return r1;
				Result r2=value(i+1,r,vX);
				if(r2.isFatalError())return r2;
				return new Result(Complex.pow(r1.val,r2.val));
			}

		// Sqrt symbol
		if(text.charAt(l)=='√'){
			interpretResult[l].submit(r,SymbolCachePair.SYMBOL_SQRT,-1);
			Result r1=value(l+1,r,vX);
			if(r1.isFatalError())return r1;
			return new Result(Complex.sqrt(r1.val));
		}

		// Brackets
		if(text.charAt(r)!=')')return new Result(1).append("Evaluator",s+" can't be calculated",l,r);
		if(text.charAt(l)=='('){
			interpretResult[l].submit(r,SymbolCachePair.SYMBOL_BRACKET,-1);
			return value(l+1,r-1,vX);
		}

		// parse function
		interpretResult[l].submit(r,SymbolCachePair.SYMBOL_FUNC,-1);
		return funcValue(l,r,vX);
	}

	private Result funcValue(int l,int r,Complex vX){ // value() for function
		String s=text.substring(l,r+1);

		// Functions
		int listPos; // the Position in Function class
		int funcID; // the ID in Function class
		int paramNum; // how many params in the function
		int leftBr; // where's the left bracket
		int exprParamNum; // how many params requires functional input

		if(funcSer[l]<0){ // not searched in list yet
			for(int i=0;i<Function.funcList.length;i++){
				if(s.startsWith(Function.funcList[i].funcName+"(")){
					//Log.i("expression","parse "+s);
					funcSer[l]=i; // found
					break;
				}
			}
		}

		listPos=funcSer[l];

		// Not found
		if(listPos<0){
			return new Result(1).append("Evaluator",s+" can't be calculated",l,r);
		}

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
				Result.setBase(Result.base);
				return new Result(new Complex(0)).append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			case Function.PREC+1:
				int prec=(int)Math.round(val[0].re);
				if(prec>Result.maxPrecision)
					return new Result(-1)
					.append("Setting","Precision too high",l,r)
					.append("Setting","The maximum precision is "+Result.maxPrecision+" digits",l,r);
				Result.precision=prec<0?Result.maxPrecision:prec;
				return new Result(new Complex(0)).append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			case Function.BASE:
				Result.setBase(10);
				return new Result(new Complex(0))
				.append("Setting","Base is set to base-"+10+" ("+ParseNumber.baseName[10]+")",l,r)
				.append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			case Function.BASE+1:
				int base=(int)Math.round(val[0].re);
				if(!(base>=2&&base<=10||base==12||base==16))
					return new Result(-1).append("Setting","Invalid base value",l,r);
				Result.setBase(base);
				return new Result(new Complex(0))
				.append("Setting","Base is set to base-"+base+" ("+ParseNumber.baseName[base]+")",l,r)
				.append("Setting","Precision is set to "+Result.precision+" digits",l,r);
			default:
		}

		return new Result(1).append("Interpreter","Function "+Function.funcList[listPos].funcName+" syntax error",l,r);
	}

	public Result value(){ // Entrance !


		isWorking=true;
		isExited=false;
		isIntegOverTolerance=false;
		isDiffOverTolerance=false;
		if(brDiff!=0){
			isExited=true;
			return new Result(1).append("Interpreter","Brackets not paired",text.length()-1,text.length()-1);
		}

		// Start working
		// 0+NaNi id a sign for "No Variable X provided" initially
		initCache();
		Result res=value(0,text.length()-1,new Complex(0,Double.NaN));
		if(isWorking&&res.val.isValid()){
			ansValue=res.val;
		}
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
		if(r1.norm2()<1E-20&&v2.norm2()<1E-20){ // 1E-10 precision
			return new Result(x2);
		}

		Complex r2=Complex.div(v2,diff(l,r,x2,r1).val); // use dir diff to speed up
		if(r2.isNaN()){
			return new Result(-1).append("Solver","Invalid initial value",l,r); // Error occurred
		}

		Complex x3;
		Result root=new Result(0);
		List<Complex> histRes=new ArrayList<>();
		double minDe=1E200;
		int minPos=-1;
		int overErrorRangeCount=0;

		for(int i=0;i<=iter;i++){ // normally no more than 20 iter., but for eq. such as x^.2, more is needed.

			if(!isWorking)return new Result(3); // Calculation Thread Halted

			Complex d1=Complex.mul(Complex.sub(x2,x1),r2);
			Complex d2=Complex.sub(r2,r1);


			x3=Complex.mul(Complex.sub(x2,Complex.div(d1,d2)),M); // Relaxation Param.
			//Log.i("expression","Solve x3="+x3);

			Complex deltaX=Complex.sub(x2,x3);
			double deltaE=deltaX.norm2();
			/*if(deltaX.norm2()<1E-20){ // 1E-10 precision
				if(v2.norm2()<1E-20){
					return new Result(x3);
				}
			}*/

			histRes.add(x3);
			if(i>0){

				if(deltaE<minDe){
					minDe=deltaE;
					minPos=i;
					overErrorRangeCount=0;
				}
				else{
					overErrorRangeCount++;
				}

				if(!x3.isFinite()||overErrorRangeCount>20){
					// return the result
					Complex res=histRes.get(minPos);
					if(minDe>1E-20||v2.norm2()>1E-18){
						root.append("Solver","Function might be ill-conditioned at certain interval",l,r);
					}
					else{
						root=new Result(res);
					}
					break;
				}

			}

			v2=value(l,r,x3).val;
			x1=x2;x2=x3;
			r1=r2;r2=Complex.div(v2,diff(l,r,x3).val); // use dir diff to speed up
			if(r2.isNaN()){ // Math Error
				if(M.re==1.0){ // 1.0 is a safely expressed double
					root=new Result(-1)
					.append("Solver","Invalid iteration occurred",l,r)
					.append("Solver","Function might be ill-conditioned at certain interval",l,r);
				}
				break;
			}
		}

		return root;
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

	// gauss's integration from x0 to x2
	// 7 nodes Gauss Quadrature
	private static final double[] gaussNodes=new double[]{ // G7 Nodes
	0.000000000000000000000000000000000e+00,
	4.058451513773971669066064120769615e-01,
	7.415311855993944398638647732807884e-01,
	9.491079123427585245261896840478513e-01
	};

	private static final double[] gaussWeights=new double[]{
	4.179591836734693877551020408163265e-01,
	3.818300505051189449503697754889751e-01,
	2.797053914892766679014677714237796e-01,
	1.294849661688696932706114326790820e-01
	};

	// 15 nodes Gauss Quadrature
	private static final double[] gaussNodes15=new double[]{ // G15 Nodes
	0.000000000000000000000000000000000e+00,
	2.011940939974345223006283033945962e-01,
	3.941513470775633698972073709810455e-01,
	5.709721726085388475372267372539106e-01,
	7.244177313601700474161860546139380e-01,
	8.482065834104272162006483207742169e-01,
	9.372733924007059043077589477102095e-01,
	9.879925180204854284895657185866126e-01
	};

	private static final double[] gaussWeights15=new double[]{
	2.025782419255612728806201999675193e-01,
	1.984314853271115764561183264438393e-01,
	1.861610000155622110268005618664228e-01,
	1.662692058169939335532008604812088e-01,
	1.395706779261543144478047945110283e-01,
	1.071592204671719350118695466858693e-01,
	7.036604748810812470926741645066734e-02,
	3.075324199611726835462839357720442e-02
	};

	Complex gaussIntegrate7(int l,int r,Complex x0,Complex x2){ // G7
		Complex lenH=new Complex(x2.re-x0.re,x2.im-x0.im);
		Complex halfH=new Complex(lenH.re/2,lenH.im/2);

		Complex t0=new Complex((x0.re+x2.re)/2,(x0.im+x2.im)/2);
		Complex t1p=new Complex(t0.re+gaussNodes[1]*halfH.re,t0.im+gaussNodes[1]*halfH.im);
		Complex t2p=new Complex(t0.re+gaussNodes[2]*halfH.re,t0.im+gaussNodes[2]*halfH.im);
		Complex t3p=new Complex(t0.re+gaussNodes[3]*halfH.re,t0.im+gaussNodes[3]*halfH.im);
		Complex t1n=new Complex(t0.re-gaussNodes[1]*halfH.re,t0.im-gaussNodes[1]*halfH.im);
		Complex t2n=new Complex(t0.re-gaussNodes[2]*halfH.re,t0.im-gaussNodes[2]*halfH.im);
		Complex t3n=new Complex(t0.re-gaussNodes[3]*halfH.re,t0.im-gaussNodes[3]*halfH.im);

		Complex y0=value(l,r,t0).val;
		Complex y1p=value(l,r,t1p).val;
		Complex y2p=value(l,r,t2p).val;
		Complex y3p=value(l,r,t3p).val;
		Complex y1n=value(l,r,t1n).val;
		Complex y2n=value(l,r,t2n).val;
		Complex y3n=value(l,r,t3n).val;

		Complex sum=new Complex(
		y0.re*gaussWeights[0]+(y1p.re+y1n.re)*gaussWeights[1]+(y2p.re+y2n.re)*gaussWeights[2]+(y3p.re+y3n.re)*gaussWeights[3],
		y0.im*gaussWeights[0]+(y1p.im+y1n.im)*gaussWeights[1]+(y2p.im+y2n.im)*gaussWeights[2]+(y3p.im+y3n.im)*gaussWeights[3]
		);
		sum=Complex.mul(sum,halfH);

		return sum;
	}

	Complex gaussIntegrate15(int l,int r,Complex x0,Complex x2){ // G15
		Complex lenH=new Complex(x2.re-x0.re,x2.im-x0.im);
		Complex halfH=new Complex(lenH.re/2,lenH.im/2);

		Complex t0=new Complex((x0.re+x2.re)/2,(x0.im+x2.im)/2);
		Complex[] tp=new Complex[7];
		Complex[] tn=new Complex[7];

		for(int i=0;i<7;i++){
			tp[i]=new Complex(t0.re+gaussNodes15[i+1]*halfH.re,t0.im+gaussNodes15[i+1]*halfH.im);
			tn[i]=new Complex(t0.re-gaussNodes15[i+1]*halfH.re,t0.im-gaussNodes15[i+1]*halfH.im);
		}

		t0=value(l,r,t0).val;
		Complex sum=new Complex(t0.re*gaussWeights15[0],t0.im*gaussWeights15[0]);
		for(int i=0;i<7;i++){
			tp[i]=value(l,r,tp[i]).val;
			tn[i]=value(l,r,tn[i]).val;
			sum=Complex.add(sum,new Complex(
				(tp[i].re+tn[i].re)*gaussWeights15[i+1],
				(tp[i].im+tn[i].im)*gaussWeights15[i+1]
			));
		}

		sum=Complex.mul(sum,halfH);

		return sum;
	}

	// integrate from x0 to x2, auto step-length
	private boolean isIntegOverTolerance=false; // only checked for once, reduce data traffic
	//private static final double logInfIntegrateStepRatio=2;
	private static final double infIntegrateStepRatio=Math.E*Math.E;
	Result adaptiveIntegrate(int l,int r,Complex x0,Complex x2,Complex lastSum,double TOL,int depth){

		//Log.i("Quadrature","Integ "+x0+" to "+x2);

		if(!isWorking)return new Result(3); // Calculation Thread Halted

		Complex x1;

		if(Double.isInfinite(x0.re)){
			return new Result(3).setVal(new Complex(0));
		}

		if(Double.isInfinite(x2.re)){
			double aRe=Math.abs(x0.re);
			double newRe;
			if(aRe<=1E5){ // Avoid precision loss due to large number
				double logNextPoint=aRe*2;
				newRe=(logNextPoint<1?Math.exp(logNextPoint):infIntegrateStepRatio*aRe);
			}
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
		sAC=gaussIntegrate15(l,r,x0,x1);
		sCB=gaussIntegrate15(l,r,x1,x2);
		sABnew=new Complex(0);
		if(sAC.isFinite())sABnew=Complex.add(sABnew,sAC);
		if(sCB.isFinite())sABnew=Complex.add(sABnew,sCB);
		abbr=Complex.sub(sAB,sABnew);

		//Log.i("Quadrature","Integ err="+abbr.norm().re);

		if(abbr.isValid()&&abbr.norm2()<200*TOL){ // Under precision limit
			return new Result(sABnew);
		}

		if(depth>=20){ // Max iteration
			Result r1=new Result(sABnew);
			if(!isIntegOverTolerance&&abbr.norm().re>1E3*TOL){
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
		if(x0.isNaN()){
			return new Result(-1).append("Quadrature","Invalid lower bound",l,r);
		}
		if(x2.isNaN()){
			return new Result(-1).append("Quadrature","Invalid upper bound",l,r);
		}
		Result check=value(l,r,x0); // Check for syntax error. Better solution ?
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

		Complex sAB=gaussIntegrate15(l,r,x0,x2);
		return adaptiveIntegrate(l,r,x0,x2,sAB,TOL*TOL,0);
	}

	// summation
	Result sum(int l,int r,Complex start,Complex end){
		double ds=start.re;
		double de=end.re;
		boolean isInfiniteSummation=(Double.isInfinite(ds)||Double.isInfinite(de));

		if(de<ds){
			return new Result(2).append("Summation","Upper bound is smaller than lower bound",l,r);
		}

		Complex sum=new Complex(0);
		Complex v=new Complex(0);
		final double TOL2=1E-16; // 8 digits expected
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

			if(isInfiniteSummation){
				if(res.val.norm2()<TOL2){
					boundCnt++;
				}else{
					boundCnt=0;
				}
				if(boundCnt>maxBoundCnt){
					break;
				}
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
				ans=Complex.mul(new Complex(n.re/m.re),ans);
				n.re-=1;
				m.re-=1;
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

	// limit (Lagrange interpolation at x0/inf)
	private static final Complex par2p=new Complex((1+Math.sqrt(2))/4);
	private static final Complex par2n=new Complex((1-Math.sqrt(2))/4);
	private static final Complex hRatio=new Complex(1+Math.sqrt(2));
	private Result limitH(int l,int r,Complex x0,Complex h){

		Complex f0;
		boolean finiteLimit=x0.isFinite();

		Complex x1,x2,x3,x4;

		if(finiteLimit){
			x1=Complex.add(x0,h);
			x2=Complex.add(x0,Complex.mul(h,hRatio));
			x3=Complex.sub(x0,h);
			x4=Complex.sub(x0,Complex.mul(h,hRatio));
		}
		else{
			double norm2=h.norm2();
			x1=new Complex(h.re/norm2,h.im/norm2);
			x2=new Complex(x1.re*2,x1.im*2);
			x3=new Complex(x1.re*3,x1.im*3);
			x4=new Complex(x1.re*4,x1.im*4);
		}

		Result r1=value(l,r,x1);if(r1.isFatalError())return r1;
		Result r2=value(l,r,x2);if(r2.isFatalError())return r2;
		Result r3=value(l,r,x3);if(r3.isFatalError())return r3;
		Result r4=value(l,r,x4);if(r4.isFatalError())return r4;

		Complex f1=r1.val;
		Complex f2=r2.val;
		Complex f3=r3.val;
		Complex f4=r4.val;

		if(finiteLimit){
			Complex f13=Complex.mul(Complex.add(f1,f3),par2p);
			Complex f24=Complex.mul(Complex.add(f2,f4),par2n);
			f0=Complex.add(f13,f24);
		}
		else{
			f0=new Complex(
			(-f1.re+24*f2.re-81*f3.re+64*f4.re)/6,
			(-f1.im+24*f2.im-81*f3.im+64*f4.im)/6
			);
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

		if(minPos<1){
			return new Result(-1).append("Evaluator","Function may be not convergent at given point",l,r);
		}
		Complex minRes=histRes.get(minPos-1);
		//Log.i("Limit","err="+minDe);

		if(minDe>1E-5){ // didn't found ?
			return new Result(-1).append("Evaluator","Function may be not convergent at given point",l,r);
		}else{ // found
			return new Result(minRes);
		}
	}

	// stop evaluation process, can be called asynchronizely
	public void stopEvaluation(){
		isWorking=false;
	}
}
