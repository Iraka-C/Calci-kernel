package complexExpression;

/*
 * Java Complex & Complex Function Class
 * Written by Iraka on 20170312
 *
 */

import java.lang.Math;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

public class Complex {
	public static Complex E=new Complex(Math.E);
	public static Complex PI=new Complex(Math.PI);
	public static Complex I=new Complex(0,1);
	public static Complex Inf=new Complex(Double.POSITIVE_INFINITY,Double.NaN);

	public double re;
	public double im;
	public Complex(double re_,double im_){
		re=re_;
		im=im_;
	}
	public Complex(double re_){
		re=re_;
		im=0;
	}
	public Complex(){
		re=Double.NaN;
		im=Double.NaN;
	}

	public static Complex add(Complex a,Complex b){
		return new Complex(a.re+b.re,a.im+b.im);
	}
	public static Complex sub(Complex a,Complex b){
		return new Complex(a.re-b.re,a.im-b.im);
	}
	public static Complex inv(Complex a){
		return new Complex(-a.re,-a.im);
	}
	public static Complex mul(Complex a,Complex b){
		return new Complex(
			a.re*b.re-a.im*b.im,
			a.re*b.im+a.im*b.re
		);
	}

	public double norm2(){
		if(Double.isInfinite(re)||Double.isInfinite(im))
			return Double.POSITIVE_INFINITY;
		return re*re+im*im;
	}

	public Complex norm(){
		/* Old implementation. There's already a hypot function

		// needs repairing: what if 1E200+1E200i ?
		// use z -> arg(z) -> cos(theta) -> norm instead.
		if(Double.isInfinite(re)||Double.isInfinite(im))
			return new Complex(Double.POSITIVE_INFINITY);
		double absRe=Math.abs(re);
		double absIm=Math.abs(im);
		if(absRe>1E150||absIm>1E150||(absRe<1E-150&&absIm<1E-150)){ // Large scales
			if(absIm>absRe){ // make absRe the bigger to increase precision
				double temp=absIm;
				absIm=absRe;
				absRe=temp;
			}
			double arg=Math.atan2(absIm,absRe);
			double cosArg=Math.cos(arg);
			return new Complex(absRe/cosArg);
		}
		return new Complex(Math.sqrt(re*re+im*im));*/
		return new Complex(Math.hypot(re,im));
	}
	public Complex arg(){
		// deal with the difference between 0 and -0
		// but better solution shall be discussed
		if(im==0){
			im=0; // -0 -> 0
			if(re==0) // the arg of 0 is not determined
				return new Complex(Double.NaN);
		}
		return new Complex(Math.atan2(im,re));
	}
	public boolean isNaN(){
		return Double.isNaN(re); // im is not cared
	}
	public boolean isValid(){ // finite complex or Complex Infinity
		return !(isDoubleFinite(re)&&Double.isNaN(im));
	}

	public static boolean isDoubleFinite(double d){
		return !(Double.isNaN(d)||Double.isInfinite(d));
	}

	public boolean isFinite(){
		return Complex.isDoubleFinite(re)&&Complex.isDoubleFinite(im);
	}

	public static Complex div(Complex a,Complex b){
		double aNorm=a.norm2();
		double bNorm=b.norm2();
		if(aNorm>0&&bNorm==0)return Inf; // pInf==nInf in complex field?
		if(Double.isInfinite(bNorm)&&Complex.isDoubleFinite(aNorm))
			return new Complex(0);
		double re=(a.re*b.re+a.im*b.im)/bNorm;
		double im=(a.im*b.re-a.re*b.im)/bNorm;
		return new Complex(re,im);
	}

	private static String significand(double oldDouble,int scale_) {
		// Original implementation of number output

		int scale=scale_>=0?scale_:0;
		String s=new BigDecimal(Double.toString(oldDouble),new MathContext(scale,RoundingMode.HALF_EVEN)).toString();
		StringBuffer sBuffer=new StringBuffer(s);
		int pos=s.indexOf('E');
		if(pos==-1)pos=s.length(); // no E found
		pos--;

		int pdot=s.indexOf('.');
		if(pdot==-1)pdot=s.length(); // no . found

		// delete redundant 0 after the dot
		while(pos>pdot&&s.charAt(pos)=='0'&&s.charAt(pos-1)>='0'&&s.charAt(pos-1)<='9'){ // you may use s, but s isn't changed
			sBuffer.deleteCharAt(pos);
			pos--;
		}
		return sBuffer.toString();
	}

	private static String doubleToString(double d){
		if(Double.isNaN(d)){
			return "nan";
		}
		if(Double.isInfinite(d)){
			return d>0?"inf":"-inf";
		}

		/*if(Result.precision<15){
			//return Complex.significand(d,Result.precision);
			return ParseNumber.toBaseString(d,2,Result.precision);
		}
		else{
			return Double.toString(d);
		}*/
		if(Result.base==10&&Result.precision==Result.maxPrecision){
			return Double.toString(d);
		}

		return ParseNumber.toBaseString(d,Result.base,Result.precision);
	}

	public String toString(){ // kind of slow !!!
		String s="";
		double threshold=(Result.precision<Result.maxPrecision?Math.pow(Result.base,-Result.precision):0);
		if(Double.isNaN(im)&&Double.isInfinite(re)){
			s=(re>0?"∞":"-∞");
		}
		else if(Math.abs(re)>threshold||Double.isNaN(re)){ // re to be shown.
			s+=doubleToString(re);

			if(isDoubleFinite(im)){
				if(Math.abs(im)>threshold){
					s+=(im>0?"+":"-");
					if(Math.abs(Math.abs(im)-1)>threshold){
						s+=doubleToString(Math.abs(im));
					}
					s+="i";
				}
			}
			else{ // inf or nan
				s+=(im<0?"":"+"); // +inf/nan -> +
				s+=doubleToString(im)+"*i";
			}
		}
		else{
			if(isDoubleFinite(im)){
				if(Math.abs(im)>threshold){
					s+=(im>0?"":"-");
					if(Math.abs(Math.abs(im)-1)>threshold){
						s+=doubleToString(Math.abs(im));
					}
					s+="i";
				}
				else{ // Nothing
					s+="0";
				}
			}
			else{ // inf nan
				s+=doubleToString(im)+"*i";
			}
		}

		return s;
	}

	//======================= Complex Functions ============================
	public static Complex exp(Complex c){
		if(c.re==Double.NEGATIVE_INFINITY)
			return new Complex(0,0);
		double norm=Math.exp(c.re);
		return new Complex(norm*Math.cos(c.im),norm*Math.sin(c.im));
	}

	public static Complex ln(Complex c){
		if(c.re==0&&c.im==0)
			return Complex.inv(Inf);
		return new Complex(Math.log(c.norm().re),c.arg().re);
	}

	public static Complex sqrt(Complex c){
		double normq=Math.sqrt(c.norm().re);
		double arg=c.arg().re;
		double sind2=Math.sqrt((1-Math.cos(arg))/2);
		double cosd2=Math.sqrt((1+Math.cos(arg))/2);
		if(arg<0)sind2=-sind2;
		return new Complex(normq*cosd2,normq*sind2);
	}

	public static Complex sin(Complex c){
		double eip=Math.exp(c.im);
		double ein=Math.exp(-c.im);
		return new Complex((eip+ein)*Math.sin(c.re)/2,(eip-ein)*Math.cos(c.re)/2);
	}

	public static Complex cos(Complex c){
		double eip=Math.exp(c.im);
		double ein=Math.exp(-c.im);
		return new Complex((eip+ein)*Math.cos(c.re)/2,(ein-eip)*Math.sin(c.re)/2);
	}

	public static Complex tan(Complex c){
		//return Complex.div(Complex.sin(c),Complex.cos(c)); // not precise enough
		double re2=c.re*2;
		double im2=c.im*2;

		double eip2=Math.exp(im2);
		double ein2=Math.exp(-im2);
		double sinhi2=(eip2-ein2)/2;
		double coshi2=(eip2+ein2)/2;

		if(Double.isInfinite(coshi2)){ // Special case
			return new Complex(0,c.im>0?1:-1);
		}

		double ratio=Math.cos(re2)+coshi2;
		double resRe=Math.sin(re2)/ratio;
		double resIm=sinhi2/ratio;
		return new Complex(resRe,resIm);
	}

	public static Complex arcsin(Complex c){
		Complex v=Complex.add(Complex.mul(c,I),Complex.sqrt(Complex.sub(new Complex(1),Complex.mul(c,c))));
		return Complex.mul(new Complex(0,-1),Complex.ln(v));
	}

	public static Complex arccos(Complex c){
		Complex v=Complex.add(c,Complex.sqrt(Complex.sub(Complex.mul(c,c),new Complex(1))));
		return Complex.mul(new Complex(0,-1),Complex.ln(v));
	}

	public static Complex arctan(Complex c){
		if(c.re==Double.POSITIVE_INFINITY)return new Complex(Math.PI/2);
		if(c.re==Double.NEGATIVE_INFINITY)return new Complex(Math.PI/2);
		/*Complex v1=Complex.ln(Complex.sub(new Complex(1),Complex.mul(c,I)));
		Complex v2=Complex.ln(Complex.add(new Complex(1),Complex.mul(c,I)));
		return Complex.mul(new Complex(0,0.5),Complex.sub(v1,v2));*/ // Old implementation

		Complex c1=new Complex(1-c.im,c.re);
		Complex c2=new Complex(1+c.im,-c.re);
		double re_=(c1.arg().re-c2.arg().re)/2;
		double im_=(Math.log(c2.norm().re)-Math.log(c1.norm().re))/2;
		return new Complex(re_,im_);
	}

	private static double[] gammaP={
		676.5203681218851,-1259.1392167224028,771.32342877765313,
		-176.61502916214059,12.507343278686905,-0.13857109526572012,
		9.9843695780195716E-6,1.5056327351493116E-7
	};
	public static Complex gamma(Complex c){ // Lanczos approximation

		if(c.re>=142.2152&&Math.abs(c.im)<=0.001){ // secure to produce complex inf
			return Complex.Inf;
		}

		Complex result;

		if(c.re<-310){ // bounded
			if(c.re==Math.floor(c.re))
				return Complex.Inf;
			else
				return new Complex(0);
		}
		else if(c.re<-0.5){ // minus complex plane
			int k=(int)Math.floor(-c.re)+1;
			result=Complex.gamma(new Complex(c.re+k,c.im));
			for(int i=k-1;i>=0;i--){ // reversed order, prevent 0/0 -> NaN
				if(!result.isFinite())break;
				result=Complex.div(result,new Complex(c.re+i,c.im));
			}
			return result;
		}
		else if(c.re<0.5){ // Reflection formula
			Complex sZ=Complex.sin(Complex.mul(Complex.PI,c));
			Complex gZ=Complex.gamma(Complex.sub(new Complex(1),c));
			result=Complex.div(Complex.PI,Complex.mul(sZ,gZ));
		}
		else{
			Complex z=new Complex(c.re-1,c.im);
			Complex x=new Complex(0.99999999999980993);

			for(int i=0;i<gammaP.length;i++){
				Complex dn=new Complex(z.re+i+1,z.im);
				x=Complex.add(x,Complex.div(new Complex(gammaP[i]),dn));
			}

			Complex t=new Complex(z.re+gammaP.length-0.5,z.im);
			result=Complex.exp(Complex.mul(new Complex(z.re+0.5,z.im),Complex.ln(t)));
			result=Complex.mul(new Complex(Math.sqrt(2*Math.PI)),result);
			result=Complex.mul(Complex.exp(Complex.inv(t)),result);
			result=Complex.mul(result,x);
		}
		return result;
	}


}
