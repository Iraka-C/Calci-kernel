package complexExpression;

/**
 * Created by Iraka Crow on 2017/3/25.
 */

class Function{

	// place a function serial here, interval 10 (accept maximum 10 params.)

	static final int EXP=10;
	static final int LN=20;
	static final int RE=30;
	static final int IM=40;
	static final int SQRT=50;
	static final int ABS=60;
	static final int NORM=70;
	static final int ARG=80;
	static final int SIN=90;
	static final int COS=100;
	static final int TAN=110;
	static final int ARCSIN=120;
	static final int ARCCOS=130;
	static final int ARCTAN=140;
	static final int GAMMA=150;
	static final int FLOOR=160;
	static final int CEIL=170;
	static final int REG=180;
	static final int CONJ=190;
	static final int RAND=200;
	static final int ROUND=210;
	static final int DIFF=220;
	static final int LIMIT=230;
	static final int EVAL=240;
	static final int FZERO=250;
	static final int INTEG=260;
	static final int SUM=270;
	static final int PERM=280;
	static final int COMB=290;
	static final int PREC=300;
	static final int COLOR=310;
	static final int BASE=320;

	static class Serial{ // function name - serial struct
		String funcName;
		int funcSerial;
		int exprParamNum; // how many param accept expression(with x) as input
		Serial(String name_,int serial_){
			funcName=name_;
			funcSerial=serial_;
			exprParamNum=0;
		}
		Serial(String name_,int serial_,int ePN_){
			funcName=name_;
			funcSerial=serial_;
			exprParamNum=ePN_;
		}
	}

	// Register name and serial pair
	// Different name may direct to same function
	static final Serial[] funcList={
	new Serial("exp",EXP),
	new Serial("ln",LN),
	new Serial("re",RE),
	new Serial("im",IM),
	new Serial("sqrt",SQRT),
	new Serial("abs",ABS),
	new Serial("norm",NORM),
	new Serial("arg",ARG),
	new Serial("sin",SIN),
	new Serial("cos",COS),
	new Serial("tan",TAN),
	new Serial("arcsin",ARCSIN),
	new Serial("arccos",ARCCOS),
	new Serial("arctan",ARCTAN),
	new Serial("gamma",GAMMA),
	new Serial("Γ",GAMMA), // greek alphabet version !
	new Serial("floor",FLOOR),
	new Serial("ceil",CEIL),
	new Serial("reg",REG),
	new Serial("conj",CONJ),
	new Serial("rand",RAND),
	new Serial("round",ROUND),
	new Serial("diff",DIFF,1),
	new Serial("limit",LIMIT,1),
	new Serial("eval",EVAL,1),
	new Serial("fzero",FZERO,1),
	new Serial("integ",INTEG,1),
	new Serial("sum",SUM,1),
	new Serial("perm",PERM),
	new Serial("comb",COMB),
	new Serial("prec",PREC),
	new Serial("color",COLOR),
	new Serial("base",BASE),
	};
}
