# Calci-kernel
A complex calculation kernel in Java (for Calci calculator)

These are complex-number classes to implement complex calculation.<br/>
The content of each class is as follows:<br/>

* Complex.java: Contains a complex number and basic complex functions
* Expression.java: Expression parser and expressional functions (such as diff, integ, fzero & so on)
* Result.java: Contains a complex result and error codes
* Function.java: Function names and serials used in Expression.java is registrated in this class
* ParseNumber.java: Number parser for different radix

***
The simplest way to parse a function and calculate its value is:
`
Expression expr = new Expression(string);
System.out.println(expr.value().val);
`
