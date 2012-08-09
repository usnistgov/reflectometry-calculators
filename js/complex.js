/*
Examples From
JavaScript: The Definitive Guide, Fourth Edition

Legal matters: these files were created by David Flanagan, and are
Copyright (c) 2001 by David Flanagan.  You may use, study, modify, and
distribute them for any purpose.  Please note that these examples are
provided "as-is" and come with no warranty of any kind.

David Flanagan
*/

/*
Modified 2012 by 
Brian Ben Maranville

additions include methods: phase, fromMagPhase, cos, sin, exp, etc.
all additions are in the public domain
*/

/*
 * Complex.js:
 * This file defines a Complex class to represent complex numbers.
 * Recall that a complex number is the sum of a real number and an
 * imaginary number, and that the imaginary number i is the
 * square root of -1.
 */

/*
 * The first step in defining a class is defining the constructor
 * function of the class. This constructor should initialize any
 * instance properties of the object. These are the essential
 * "state variables" that make each instance of the class different.
 */
function Complex(real, imaginary) {
    this.x = real;       // The real part of the number
    this.y = imaginary;  // The imaginary part of the number
}

/*
 * The second step in defining a class is defining its instance
 * methods (and possibly other properties) in the prototype object
 * of the constructor. Any properties defined in this object will
 * be inherited by all instances of the class. Note that instance
 * methods operate implicitly on the this keyword. For many methods,
 * no other arguments are needed.
 */

// Return the magnitude of a complex number. This is defined
// as its distance from the origin (0,0) of the complex plane.
Complex.prototype.magnitude = function() {
    return Math.sqrt(this.x*this.x + this.y*this.y);
};

// Return the phase of a complex number. This is defined
// as its angle (counterclockwise) from the real axis of the complex plane.
Complex.prototype.phase = function() {
    return Math.atan2(this.y, this.x);
};
// Return a complex number that is the negative of this one.
Complex.prototype.negative = function() {
    return new Complex(-this.x, -this.y);
};

// Return an exact copy of the complex number
Complex.prototype.copy = function() {
    return new Complex(this.x, this.y);
};

Complex.prototype.inverse = function() {
    var denom = Math.pow(this.x, 2) + Math.pow(this.y, 2);
    return new Complex( this.x / denom, -this.y / denom );
};

Complex.prototype.conjugate = function() {
    return new Complex( this.x, -this.y );
};

// Convert a Complex object to a string in a useful way.
// This is invoked when a Complex object is used as a string.
Complex.prototype.toString = function() {
    return "{" + this.x + "," + this.y + "}";
};

// Return the real portion of a complex number. This function
// is invoked when a Complex object is treated as a primitive value.
Complex.prototype.valueOf = function() { return this.x; }

/*
 * The third step in defining a class is to define class methods,
 * constants, and any needed class properties as properties of the
 * constructor function itself (instead of as properties of the
 * prototype object of the constructor). Note that class methods
 * do not use the this keyword: they operate only on their arguments.
 */

// Add two complex numbers and return the result.
Complex.add = function (a, b) {
    if (!(a instanceof Complex)) a = new Complex(a, 0);
    if (!(b instanceof Complex)) b = new Complex(b, 0);
    return new Complex(a.x + b.x, a.y + b.y);
};

// Subtract one complex number from another.
Complex.subtract = function (a, b) {
    if (!(a instanceof Complex)) a = new Complex(a, 0);
    if (!(b instanceof Complex)) b = new Complex(b, 0);
    return new Complex(a.x - b.x, a.y - b.y);
};

// Multiply two complex numbers and return the product.
Complex.multiply = function(a, b) {
    if (!(a instanceof Complex)) a = new Complex(a, 0);
    if (!(b instanceof Complex)) b = new Complex(b, 0);
    return new Complex(a.x * b.x - a.y * b.y,
                       a.x * b.y + a.y * b.x);
};

// Here are some useful predefined complex numbers.
// They are defined as class properties, where they can be used as
// "constants." (Note, though, that they are not actually read-only.)
Complex.zero = new Complex(0,0);
Complex.one = new Complex(1,0);
Complex.i = new Complex(0,1);
Complex.half = new Complex(0.5, 0);

Complex.exp = function(a) {
    var re = Math.exp(a.x);
    return new Complex(re * Math.cos(a.y), re * Math.sin(a.y));
};

Complex.fromMagPhase = function(mag, phase) {
    var x = Math.cos(phase);
    var y = Math.sin(phase);
    return new Complex( mag * Math.cos(phase), mag * Math.sin(phase) );
};

Complex.cos = function(a) {
    var t1 = Complex.exp(Complex.multiply(a, Complex.i));
    var t2 = Complex.exp(Complex.multiply(a, Complex.i.negative()));
    var d = new Complex(2, 0);
    return Complex.multiply(Complex.add(t1, t2), d.inverse());
    // (e^ia + e^-ia)/2
};

Complex.sin = function(a) {
    var t1 = Complex.exp(Complex.multiply(a, Complex.i));
    var t2 = Complex.exp(Complex.multiply(a, Complex.i.negative()));
    var d = new Complex(0, 2);
    return Complex.multiply(Complex.subtract(t1, t2), d.inverse());
    // (e^ia - e^-ia)/2i
};

Complex.sqrt = function(a) {
    if (!(a instanceof Complex)) a = new Complex(a, 0);
    var mag = Math.sqrt(a.magnitude());
    var phase = 0.5 * a.phase();
    return new Complex( mag * Math.cos(phase), mag * Math.sin(phase));
    // sqrt(r*e^ia) = sqrt(r) * e^i(a/2)
};

Complex.sum = function(alist) {
    var result = new Complex(0,0);
    for (var i = 0; i<alist.length; i++) {
        result = Complex.add(result, alist[i]);
    }
    return result
};

Complex.cosh = function(a) {
    return Complex.multiply(Complex.half, Complex.add( Complex.exp(a), Complex.exp(a.negative()) ) );
    // cosh(x) = 0.5*(e^x + e^-x)
};

Complex.sinh = function(a) {
    return Complex.multiply(Complex.half, Complex.subtract( Complex.exp(a), Complex.exp(a.negative()) ) );
    // sinh(x) = 0.5 * (e^x - e^-x);
};
    
    
