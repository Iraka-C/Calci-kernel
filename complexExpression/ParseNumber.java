import java.math.BigDecimal;

/**
 * Created by Iraka Crow on 2017/4/12.
 */

public class ParseNumber{

	private static final String baseSymbol="₀₁₂₃₄₅₆₇₈₉₉₉ₑₑₑₑₓ";

	// is this character a base notation?
	public static boolean isBaseSymbol(char c){
		return baseSymbol.indexOf(c)>=0;
	}

	// get a number digit from char under certain base
	private static int getDigit(char c,int base) throws NumberFormatException{
		int digit;
		if(c>='0'&&c<='9'){ // in base 10
			digit=c-'0';
		}
		else if(c>='A'&&c<='F'){ // in base 16
			digit=c-'A'+10;
		}
		else{ // not a valid digit
			throw new NumberFormatException();
		}
		if(digit>=base)throw new NumberFormatException();
		return digit;
	}

	private static double parseRaw(String s,int base)throws NumberFormatException{
		int dotPos;
		for(dotPos=0;dotPos<s.length();dotPos++){
			if(s.charAt(dotPos)=='.')
				break;
		}

		double frac=0;
		double digitBase=1;
		for(int i=dotPos-1;i>=0;i--){
			frac+=getDigit(s.charAt(i),base)*digitBase;
			digitBase*=base;
		}
		digitBase=1./base;
		for(int i=dotPos+1;i<s.length();i++){
			frac+=getDigit(s.charAt(i),base)*digitBase;
			digitBase/=base;
		}

		return frac;
	}

	// parse a float number presentation under certain base
	public static double parse(String s)throws NumberFormatException{

		int base=0;
		int baseSymbolPos=-1;
		for(int i=0;i<s.length();i++){
			base=baseSymbol.indexOf(s.charAt(i));
			if(base>0){ // base 0 if invalid
				baseSymbolPos=i;
				break;
			}
		}
		if(baseSymbolPos==0)throw new NumberFormatException();
		if(baseSymbolPos<0)return parseCompat(s); // parse the string as a compat format


		int exp;
		if(baseSymbolPos==s.length()-1){
			exp=0;
		}
		else{
			exp=Integer.parseInt(s.substring(baseSymbolPos+1));
		}
		double frac=parseRaw(s.substring(0,baseSymbolPos),base);

		return frac*Math.pow(base,exp);
	}

	public static double parseCompat(String s)throws NumberFormatException{
		int baseDivSymbolPos=s.indexOf('~');
		if(baseDivSymbolPos<=0||baseDivSymbolPos>=s.length()-1)
			throw new NumberFormatException();

		int base=Integer.parseInt(s.substring(baseDivSymbolPos+1));
		if(!(base>=2&&base<=10||base==12||base==16)) // base not supported
			throw new NumberFormatException();

		return parseRaw(s.substring(0,baseDivSymbolPos),base);
	}

	// without scientific display mode
	private static final String numSymbol="0123456789ABCDEF";
	private static String toPositiveRawBaseString(double d,long base,int prec){ // d_>0
		int[] digits=new int[100];
		int intDigitNum=(int)Math.floor(Math.log(d)/Math.log(base))+1;
		if(intDigitNum<0)intDigitNum=0;

		long intPart=(long)Math.floor(d);
		double fracPart=d-intPart;

		for(int i=intDigitNum;i>=0;i--){
			digits[i]=(int)(intPart%base);
			intPart/=base;
		}

		for(int i=intDigitNum+1;i<=prec+1;i++){
			fracPart*=base;
			digits[i]=(int)Math.floor(fracPart);
			fracPart-=digits[i];
		}

		if(digits[prec+1]*2>=base){ // carry-over
			digits[prec]++;
			for(int i=prec;i>0;i--){
				if(digits[i]==base){
					digits[i]=0;
					digits[i-1]++;
				}
				else{
					break;
				}
			}
		}

		int maxNonZeroPos; // omit 0 digits after dot & at the rear of a notation
		for(maxNonZeroPos=prec;maxNonZeroPos>=0;maxNonZeroPos--){
			if(maxNonZeroPos<=intDigitNum||digits[maxNonZeroPos]>0)
				break;
		}

		String result="";
		for(int i=0;i<=maxNonZeroPos;i++){ // construct string
			if(!(intDigitNum>0&&i==0&&digits[0]==0))
				result+=numSymbol.charAt(digits[i]);
			if(i==intDigitNum&&i<maxNonZeroPos)
				result+='.';
		}

		for(int i=maxNonZeroPos+1;i<=intDigitNum;i++)
			result+='0';

		return result;
	}

	// print the value of d under base to a string with prec digits' precision
	public static String toBaseString(double d_,int base,int prec){
		if(Double.isNaN(d_))return "nan";
		if(d_==Double.POSITIVE_INFINITY)return "inf";
		if(d_==Double.NEGATIVE_INFINITY)return "-inf";

		String negativeSymbol=(d_>=0?"":"-");
		double d=Math.abs(d_);
		double maxPreciseValue=Math.pow(base,prec);
		double minPreciseValue=Math.pow(base,-prec);

		String result;

		if(d<maxPreciseValue&&d>minPreciseValue){ // able to express under fixed precision
			result=negativeSymbol+toPositiveRawBaseString(d,base,prec);
			if(base!=10){
				result+=baseSymbol.charAt(base);
			}
		}
		else{ // need scientific notation
			double fracPart=d;
			int digitExp=0;
			while(fracPart>=base){
				digitExp++;
				fracPart/=base;
			}
			while(fracPart<1){
				digitExp--;
				fracPart*=base;
			}

			if(base==10){
				result=negativeSymbol+toPositiveRawBaseString(fracPart,base,prec);
				result+="E";
			}
			else{
				result=negativeSymbol+toPositiveRawBaseString(fracPart,base,prec);
				result+=baseSymbol.charAt(base);
			}

			result+=digitExp;
		}
		return result;
	}

	public static final String[] baseName=new String[]{
	"---","---","Binary","Ternary",
	"Quaternary","Quinary","Senary","Septenary",
	"Octal","Nonary","Decimal","Undecimal",
	"Duodecimal","Tridecimal","Tetradecimal","Pentadecimal",
	"Hexadecimal"
	};
}
