// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "legendre.h"

int legendre_order_1d[] = {

0,
1,
2,
3,
4,
5,
6,
7,
8,
9,
10,
11,
12,
13,
14,
15,
16,
17,
18,
19,
20,
21,
22,
23,
24,
25,
26,
27,
28,
29,
};


static double legendre_fn_0(double _x) {
    long double x = _x;
    return pow(2,(1.0/2.0))/2;
}

static double legendre_fn_1(double _x) {
    long double x = _x;
    return x*pow(6,(1.0/2.0))/2;
}

static double legendre_fn_2(double _x) {
    long double x = _x;
    return -pow(10,(1.0/2.0))/4 + 3*pow(10,(1.0/2.0))*((x)*(x))/4;
}

static double legendre_fn_3(double _x) {
    long double x = _x;
    return -3*x*pow(14,(1.0/2.0))*(1 - 5*((x)*(x))/3)/4;
}

static double legendre_fn_4(double _x) {
    long double x = _x;
    return -45*pow(2,(1.0/2.0))*((x)*(x))*(1 - 7*((x)*(x))/6)/8 + 9*pow(2,(1.0/2.0))/16;
}

static double legendre_fn_5(double _x) {
    long double x = _x;
    return 15*x*pow(22,(1.0/2.0))*(1 - 14*((x)*(x))*(1 - 9*((x)*(x))/10)/3)/16;
}

static double legendre_fn_6(double _x) {
    long double x = _x;
    return 105*pow(26,(1.0/2.0))*((x)*(x))*(1 - 3*((x)*(x))*(1 - 11*((x)*(x))/15))/32 - 5*pow(26,(1.0/2.0))/32;
}

static double legendre_fn_7(double _x) {
    long double x = _x;
    return -35*x*pow(30,(1.0/2.0))*(1 - 9*((x)*(x))*(1 - 11*((x)*(x))*(1 - 13*((x)*(x))/21)/5))/32;
}

static double legendre_fn_8(double _x) {
    long double x = _x;
    return -315*pow(34,(1.0/2.0))*((x)*(x))*(1 - 11*((x)*(x))*(1 - 26*((x)*(x))*(1 - 15*((x)*(x))/28)/15)/2)/64 + 35*pow(34,(1.0/2.0))/256;
}

static double legendre_fn_9(double _x) {
    long double x = _x;
    return 315*x*pow(38,(1.0/2.0))*(1 - 44*((x)*(x))*(1 - 39*((x)*(x))*(1 - 10*((x)*(x))*(1 - 17*((x)*(x))/36)/7)/10)/3)/256;
}

static double legendre_fn_10(double _x) {
    long double x = _x;
    return 3465*pow(42,(1.0/2.0))*((x)*(x))*(1 - 26*((x)*(x))*(1 - 3*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))/45)/14))/3)/512 - 63*pow(42,(1.0/2.0))/512;
}

static double legendre_fn_11(double _x) {
    long double x = _x;
    return -693*x*pow(46,(1.0/2.0))*(1 - 65*((x)*(x))*(1 - 6*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))*(1 - 21*((x)*(x))/55)/18)/7))/3)/512;
}

static double legendre_fn_12(double _x) {
    long double x = _x;
    return -45045*pow(2,(1.0/2.0))*((x)*(x))*(1 - 25*((x)*(x))*(1 - 68*((x)*(x))*(1 - 57*((x)*(x))*(1 - 14*((x)*(x))*(1 - 23*((x)*(x))/66)/15)/28)/15)/2)/1024 + 1155*pow(2,(1.0/2.0))/2048;
}

static double legendre_fn_13(double _x) {
    long double x = _x;
    return 9009*x*pow(6,(1.0/2.0))*(1 - 30*((x)*(x))*(1 - 17*((x)*(x))*(1 - 76*((x)*(x))*(1 - 7*((x)*(x))*(1 - 46*((x)*(x))*(1 - 25*((x)*(x))/78)/55)/4)/21)/2))/2048;
}

static double legendre_fn_14(double _x) {
    long double x = _x;
    return 45045*pow(58,(1.0/2.0))*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))*(1 - 3*((x)*(x))*(1 - 23*((x)*(x))*(1 - 25*((x)*(x))*(1 - 27*((x)*(x))/91)/33)/15))/3))/4096 - 429*pow(58,(1.0/2.0))/4096;
}

static double legendre_fn_15(double _x) {
    long double x = _x;
    return -6435*x*pow(62,(1.0/2.0))*(1 - 119*((x)*(x))*(1 - 57*((x)*(x))*(1 - 5*((x)*(x))*(1 - 23*((x)*(x))*(1 - 15*((x)*(x))*(1 - 9*((x)*(x))*(1 - 29*((x)*(x))/105)/13)/11)/9))/5)/3)/4096;
}

static double legendre_fn_16(double _x) {
    long double x = _x;
    return -109395*pow(66,(1.0/2.0))*((x)*(x))*(1 - 133*((x)*(x))*(1 - 42*((x)*(x))*(1 - 115*((x)*(x))*(1 - 20*((x)*(x))*(1 - 27*((x)*(x))*(1 - 58*((x)*(x))*(1 - 31*((x)*(x))/120)/91)/22)/9)/28)/5)/6)/8192 + 6435*pow(66,(1.0/2.0))/65536;
}

static double legendre_fn_17(double _x) {
    long double x = _x;
    return 109395*x*pow(70,(1.0/2.0))*(1 - 152*((x)*(x))*(1 - 147*((x)*(x))*(1 - 46*((x)*(x))*(1 - 125*((x)*(x))*(1 - 108*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 33*((x)*(x))/136)/105)/26)/55)/36)/7)/10)/3)/65536;
}

static double legendre_fn_18(double _x) {
    long double x = _x;
    return 2078505*pow(74,(1.0/2.0))*((x)*(x))*(1 - 28*((x)*(x))*(1 - 161*((x)*(x))*(1 - 75*((x)*(x))*(1 - 3*((x)*(x))*(1 - 58*((x)*(x))*(1 - 93*((x)*(x))*(1 - 11*((x)*(x))*(1 - 35*((x)*(x))/153)/20)/91)/33))/14)/15))/131072 - 12155*pow(74,(1.0/2.0))/131072;
}

static double legendre_fn_19(double _x) {
    long double x = _x;
    return -230945*x*pow(78,(1.0/2.0))*(1 - 63*((x)*(x))*(1 - 92*((x)*(x))*(1 - 25*((x)*(x))*(1 - 9*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 33*((x)*(x))*(1 - 35*((x)*(x))*(1 - 37*((x)*(x))/171)/68)/35)/39)/11)/2)/3)/5))/131072;
}

static double legendre_fn_20(double _x) {
    long double x = _x;
    return -4849845*pow(82,(1.0/2.0))*((x)*(x))*(1 - 69*((x)*(x))*(1 - 40*((x)*(x))*(1 - 27*((x)*(x))*(1 - 58*((x)*(x))*(1 - 155*((x)*(x))*(1 - 132*((x)*(x))*(1 - 7*((x)*(x))*(1 - 74*((x)*(x))*(1 - 39*((x)*(x))/190)/153)/8)/91)/66)/15)/4)/3)/2)/262144 + 46189*pow(82,(1.0/2.0))/524288;
}

static double legendre_fn_21(double _x) {
    long double x = _x;
    return 969969*x*pow(86,(1.0/2.0))*(1 - 230*((x)*(x))*(1 - 45*((x)*(x))*(1 - 72*((x)*(x))*(1 - 203*((x)*(x))*(1 - 186*((x)*(x))*(1 - 55*((x)*(x))*(1 - 4*((x)*(x))*(1 - 111*((x)*(x))*(1 - 26*((x)*(x))*(1 - 41*((x)*(x))/210)/57)/136)/3)/26)/55)/36)/7)/2)/3)/524288;
}

static double legendre_fn_22(double _x) {
    long double x = _x;
    return 66927861*pow(10,(1.0/2.0))*((x)*(x))*(1 - 125*((x)*(x))*(1 - 81*((x)*(x))*(1 - 58*((x)*(x))*(1 - 217*((x)*(x))*(1 - 3*((x)*(x))*(1 - 25*((x)*(x))*(1 - 37*((x)*(x))*(1 - 13*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))/231)/95)/17)/30)/13))/45)/7)/5)/3)/1048576 - 264537*pow(10,(1.0/2.0))/1048576;
}

static double legendre_fn_23(double _x) {
    long double x = _x;
    return -2028117*x*pow(94,(1.0/2.0))*(1 - 275*((x)*(x))*(1 - 27*((x)*(x))*(1 - 87*((x)*(x))*(1 - 62*((x)*(x))*(1 - 21*((x)*(x))*(1 - 35*((x)*(x))*(1 - 37*((x)*(x))*(1 - 39*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))*(1 - 45*((x)*(x))/253)/105)/57)/34)/21)/13)/5)/9)/7))/3)/1048576;
}

static double legendre_fn_24(double _x) {
    long double x = _x;
    return -354920475*pow(2,(1.0/2.0))*((x)*(x))*(1 - 99*((x)*(x))*(1 - 58*((x)*(x))*(1 - 279*((x)*(x))*(1 - 88*((x)*(x))*(1 - 245*((x)*(x))*(1 - 222*((x)*(x))*(1 - 13*((x)*(x))*(1 - 164*((x)*(x))*(1 - 129*((x)*(x))*(1 - 30*((x)*(x))*(1 - 47*((x)*(x))/276)/77)/190)/153)/8)/91)/66)/15)/28)/3)/2)/2097152 + 4732273*pow(2,(1.0/2.0))/8388608;
}

static double legendre_fn_25(double _x) {
    long double x = _x;
    return 16900975*x*pow(102,(1.0/2.0))*(1 - 108*((x)*(x))*(1 - 319*((x)*(x))*(1 - 310*((x)*(x))*(1 - 33*((x)*(x))*(1 - 56*((x)*(x))*(1 - 259*((x)*(x))*(1 - 78*((x)*(x))*(1 - 205*((x)*(x))*(1 - 172*((x)*(x))*(1 - 9*((x)*(x))*(1 - 94*((x)*(x))*(1 - 49*((x)*(x))/300)/253)/14)/171)/136)/35)/78)/11)/4)/21)/10))/8388608;
}

static double legendre_fn_26(double _x) {
    long double x = _x;
    return 456326325*pow(106,(1.0/2.0))*((x)*(x))*(1 - 58*((x)*(x))*(1 - 341*((x)*(x))*(1 - 165*((x)*(x))*(1 - 7*((x)*(x))*(1 - 148*((x)*(x))*(1 - 3*((x)*(x))*(1 - 41*((x)*(x))*(1 - 215*((x)*(x))*(1 - 18*((x)*(x))*(1 - 47*((x)*(x))*(1 - 49*((x)*(x))*(1 - 51*((x)*(x))/325)/138)/77)/19)/153)/20))/33))/14)/15))/16777216 - 1300075*pow(106,(1.0/2.0))/16777216;
}

static double legendre_fn_27(double _x) {
    long double x = _x;
    return -35102025*x*pow(110,(1.0/2.0))*(1 - 377*((x)*(x))*(1 - 186*((x)*(x))*(1 - 121*((x)*(x))*(1 - 175*((x)*(x))*(1 - 333*((x)*(x))*(1 - 4*((x)*(x))*(1 - 41*((x)*(x))*(1 - 129*((x)*(x))*(1 - 25*((x)*(x))*(1 - 94*((x)*(x))*(1 - 147*((x)*(x))*(1 - 17*((x)*(x))*(1 - 53*((x)*(x))/351)/50)/253)/105)/19)/68)/15))/55)/18)/7)/5)/3)/16777216;
}

static double legendre_fn_28(double _x) {
    long double x = _x;
    return -1017958725*pow(114,(1.0/2.0))*((x)*(x))*(1 - 403*((x)*(x))*(1 - 132*((x)*(x))*(1 - 55*((x)*(x))*(1 - 74*((x)*(x))*(1 - 117*((x)*(x))*(1 - 328*((x)*(x))*(1 - 301*((x)*(x))*(1 - 30*((x)*(x))*(1 - 47*((x)*(x))*(1 - 28*((x)*(x))*(1 - 51*((x)*(x))*(1 - 106*((x)*(x))*(1 - 55*((x)*(x))/378)/325)/92)/33)/38)/17)/120)/91)/22)/9)/4)/5)/6)/33554432 + 5014575*pow(114,(1.0/2.0))/67108864;
}

static double legendre_fn_29(double _x) {
    long double x = _x;
    return 145422675*x*pow(118,(1.0/2.0))*(1 - 434*((x)*(x))*(1 - 429*((x)*(x))*(1 - 20*((x)*(x))*(1 - 407*((x)*(x))*(1 - 78*((x)*(x))*(1 - 123*((x)*(x))*(1 - 344*((x)*(x))*(1 - 315*((x)*(x))*(1 - 94*((x)*(x))*(1 - 7*((x)*(x))*(1 - 204*((x)*(x))*(1 - 53*((x)*(x))*(1 - 110*((x)*(x))*(1 - 57*((x)*(x))/406)/351)/100)/253)/6)/57)/136)/105)/26)/11)/36))/10)/3)/67108864;
}


shape_fn_t legendre_fn_tab_1d[] = {

legendre_fn_0,
legendre_fn_1,
legendre_fn_2,
legendre_fn_3,
legendre_fn_4,
legendre_fn_5,
legendre_fn_6,
legendre_fn_7,
legendre_fn_8,
legendre_fn_9,
legendre_fn_10,
legendre_fn_11,
legendre_fn_12,
legendre_fn_13,
legendre_fn_14,
legendre_fn_15,
legendre_fn_16,
legendre_fn_17,
legendre_fn_18,
legendre_fn_19,
legendre_fn_20,
legendre_fn_21,
legendre_fn_22,
legendre_fn_23,
legendre_fn_24,
legendre_fn_25,
legendre_fn_26,
legendre_fn_27,
legendre_fn_28,
legendre_fn_29,
};



static double legendre_der_0(double _x) {
    long double x = _x;
    return 0;
}

static double legendre_der_1(double _x) {
    long double x = _x;
    return pow(6,(1.0/2.0))/2;
}

static double legendre_der_2(double _x) {
    long double x = _x;
    return 3*x*pow(10,(1.0/2.0))/2;
}

static double legendre_der_3(double _x) {
    long double x = _x;
    return -3*pow(14,(1.0/2.0))/4 + 15*pow(14,(1.0/2.0))*((x)*(x))/4;
}

static double legendre_der_4(double _x) {
    long double x = _x;
    return -45*x*pow(2,(1.0/2.0))*(1 - 7*((x)*(x))/3)/4;
}

static double legendre_der_5(double _x) {
    long double x = _x;
    return -105*pow(22,(1.0/2.0))*((x)*(x))*(1 - 3*((x)*(x))/2)/8 + 15*pow(22,(1.0/2.0))/16;
}

static double legendre_der_6(double _x) {
    long double x = _x;
    return 105*x*pow(26,(1.0/2.0))*(1 - 6*((x)*(x))*(1 - 11*((x)*(x))/10))/16;
}

static double legendre_der_7(double _x) {
    long double x = _x;
    return 945*pow(30,(1.0/2.0))*((x)*(x))*(1 - 11*((x)*(x))*(1 - 13*((x)*(x))/15)/3)/32 - 35*pow(30,(1.0/2.0))/32;
}

static double legendre_der_8(double _x) {
    long double x = _x;
    return -315*x*pow(34,(1.0/2.0))*(1 - 11*((x)*(x))*(1 - 13*((x)*(x))*(1 - 5*((x)*(x))/7)/5))/32;
}

static double legendre_der_9(double _x) {
    long double x = _x;
    return -3465*pow(38,(1.0/2.0))*((x)*(x))*(1 - 13*((x)*(x))*(1 - 2*((x)*(x))*(1 - 17*((x)*(x))/28))/2)/64 + 315*pow(38,(1.0/2.0))/256;
}

static double legendre_der_10(double _x) {
    long double x = _x;
    return 3465*x*pow(42,(1.0/2.0))*(1 - 52*((x)*(x))*(1 - 9*((x)*(x))*(1 - 34*((x)*(x))*(1 - 19*((x)*(x))/36)/21)/2)/3)/256;
}

static double legendre_der_11(double _x) {
    long double x = _x;
    return 45045*pow(46,(1.0/2.0))*((x)*(x))*(1 - 10*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))*(1 - 7*((x)*(x))/15)/14)/5))/512 - 693*pow(46,(1.0/2.0))/512;
}

static double legendre_der_12(double _x) {
    long double x = _x;
    return -45045*x*pow(2,(1.0/2.0))*(1 - 25*((x)*(x))*(1 - 34*((x)*(x))*(1 - 19*((x)*(x))*(1 - 7*((x)*(x))*(1 - 23*((x)*(x))/55)/6)/7)/5))/512;
}

static double legendre_der_13(double _x) {
    long double x = _x;
    return -405405*pow(6,(1.0/2.0))*((x)*(x))*(1 - 85*((x)*(x))*(1 - 76*((x)*(x))*(1 - 9*((x)*(x))*(1 - 46*((x)*(x))*(1 - 25*((x)*(x))/66)/45)/4)/15)/6)/1024 + 9009*pow(6,(1.0/2.0))/2048;
}

static double legendre_der_14(double _x) {
    long double x = _x;
    return 45045*x*pow(58,(1.0/2.0))*(1 - 34*((x)*(x))*(1 - 19*((x)*(x))*(1 - 4*((x)*(x))*(1 - 23*((x)*(x))*(1 - 10*((x)*(x))*(1 - 9*((x)*(x))/26)/11)/12))/2))/2048;
}

static double legendre_der_15(double _x) {
    long double x = _x;
    return 765765*pow(62,(1.0/2.0))*((x)*(x))*(1 - 19*((x)*(x))*(1 - 7*((x)*(x))*(1 - 23*((x)*(x))*(1 - 5*((x)*(x))*(1 - 9*((x)*(x))*(1 - 29*((x)*(x))/91)/11)/3)/7)))/4096 - 6435*pow(62,(1.0/2.0))/4096;
}

static double legendre_der_16(double _x) {
    long double x = _x;
    return -109395*x*pow(66,(1.0/2.0))*(1 - 133*((x)*(x))*(1 - 63*((x)*(x))*(1 - 115*((x)*(x))*(1 - 25*((x)*(x))*(1 - 81*((x)*(x))*(1 - 29*((x)*(x))*(1 - 31*((x)*(x))/105)/39)/55)/9)/21)/5)/3)/4096;
}

static double legendre_der_17(double _x) {
    long double x = _x;
    return -2078505*pow(70,(1.0/2.0))*((x)*(x))*(1 - 49*((x)*(x))*(1 - 46*((x)*(x))*(1 - 125*((x)*(x))*(1 - 12*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 11*((x)*(x))/40)/91)/22)/5)/28)/5)/2)/8192 + 109395*pow(70,(1.0/2.0))/65536;
}

static double legendre_der_18(double _x) {
    long double x = _x;
    return 2078505*x*pow(74,(1.0/2.0))*(1 - 56*((x)*(x))*(1 - 161*((x)*(x))*(1 - 50*((x)*(x))*(1 - 15*((x)*(x))*(1 - 116*((x)*(x))*(1 - 31*((x)*(x))*(1 - 22*((x)*(x))*(1 - 35*((x)*(x))/136)/35)/26)/55)/4)/7)/10))/65536;
}

static double legendre_der_19(double _x) {
    long double x = _x;
    return 43648605*pow(78,(1.0/2.0))*((x)*(x))*(1 - 92*((x)*(x))*(1 - 35*((x)*(x))*(1 - 81*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 99*((x)*(x))*(1 - 7*((x)*(x))*(1 - 37*((x)*(x))/153)/12)/91)/33)/9)/14)/3)/3)/131072 - 230945*pow(78,(1.0/2.0))/131072;
}

static double legendre_der_20(double _x) {
    long double x = _x;
    return -4849845*x*pow(82,(1.0/2.0))*(1 - 69*((x)*(x))*(1 - 20*((x)*(x))*(1 - 9*((x)*(x))*(1 - 29*((x)*(x))*(1 - 31*((x)*(x))*(1 - 22*((x)*(x))*(1 - ((x)*(x))*(1 - 37*((x)*(x))*(1 - 13*((x)*(x))/57)/68))/13)/11)/6))))/131072;
}

static double legendre_der_21(double _x) {
    long double x = _x;
    return -111546435*pow(86,(1.0/2.0))*((x)*(x))*(1 - 75*((x)*(x))*(1 - 72*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 5*((x)*(x))*(1 - 20*((x)*(x))*(1 - 37*((x)*(x))*(1 - 26*((x)*(x))*(1 - 41*((x)*(x))/190)/51)/40)/13)/2)/15)/4)/5)/2)/262144 + 969969*pow(86,(1.0/2.0))/524288;
}

static double legendre_der_22(double _x) {
    long double x = _x;
    return 66927861*x*pow(10,(1.0/2.0))*(1 - 250*((x)*(x))*(1 - 243*((x)*(x))*(1 - 232*((x)*(x))*(1 - 217*((x)*(x))*(1 - 18*((x)*(x))*(1 - 175*((x)*(x))*(1 - 148*((x)*(x))*(1 - 117*((x)*(x))*(1 - 82*((x)*(x))*(1 - 43*((x)*(x))/210)/171)/136)/105)/78)/5)/36)/21)/10)/3)/524288;
}

static double legendre_der_23(double _x) {
    long double x = _x;
    return 557732175*pow(94,(1.0/2.0))*((x)*(x))*(1 - 45*((x)*(x))*(1 - 87*((x)*(x))*(1 - 62*((x)*(x))*(1 - 77*((x)*(x))*(1 - 35*((x)*(x))*(1 - 185*((x)*(x))*(1 - 13*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))*(1 - 15*((x)*(x))/77)/95)/51)/10)/91)/11)/15)/7)/5))/1048576 - 2028117*pow(94,(1.0/2.0))/1048576;
}

static double legendre_der_24(double _x) {
    long double x = _x;
    return -354920475*x*pow(2,(1.0/2.0))*(1 - 99*((x)*(x))*(1 - 29*((x)*(x))*(1 - 93*((x)*(x))*(1 - 22*((x)*(x))*(1 - 49*((x)*(x))*(1 - 37*((x)*(x))*(1 - 13*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))*(1 - 3*((x)*(x))*(1 - 47*((x)*(x))/253)/7)/57)/34)/7)/13)/11)/3)/7)))/1048576;
}

static double legendre_der_25(double _x) {
    long double x = _x;
    return -1368978975*pow(102,(1.0/2.0))*((x)*(x))*(1 - 319*((x)*(x))*(1 - 62*((x)*(x))*(1 - 297*((x)*(x))*(1 - 56*((x)*(x))*(1 - 259*((x)*(x))*(1 - 18*((x)*(x))*(1 - 41*((x)*(x))*(1 - 172*((x)*(x))*(1 - 27*((x)*(x))*(1 - 94*((x)*(x))*(1 - 49*((x)*(x))/276)/231)/38)/153)/24)/7)/66)/9)/28)/3)/6)/2097152 + 16900975*pow(102,(1.0/2.0))/8388608;
}

static double legendre_der_26(double _x) {
    long double x = _x;
    return 456326325*x*pow(106,(1.0/2.0))*(1 - 116*((x)*(x))*(1 - 341*((x)*(x))*(1 - 110*((x)*(x))*(1 - 35*((x)*(x))*(1 - 296*((x)*(x))*(1 - 7*((x)*(x))*(1 - 82*((x)*(x))*(1 - 215*((x)*(x))*(1 - 20*((x)*(x))*(1 - 47*((x)*(x))*(1 - 98*((x)*(x))*(1 - 17*((x)*(x))/100)/253)/70)/19)/136)/35)/2)/55)/4)/7)/10))/8388608;
}

static double legendre_der_27(double _x) {
    long double x = _x;
    return 13233463425*pow(110,(1.0/2.0))*((x)*(x))*(1 - 62*((x)*(x))*(1 - 121*((x)*(x))*(1 - 25*((x)*(x))*(1 - 37*((x)*(x))*(1 - 52*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))*(1 - 25*((x)*(x))*(1 - 94*((x)*(x))*(1 - 7*((x)*(x))*(1 - 17*((x)*(x))*(1 - 53*((x)*(x))/325)/46)/11)/95)/17)/20)/13)/11)/5)/2)/5))/16777216 - 35102025*pow(110,(1.0/2.0))/16777216;
}

static double legendre_der_28(double _x) {
    long double x = _x;
    return -1017958725*x*pow(114,(1.0/2.0))*(1 - 403*((x)*(x))*(1 - 198*((x)*(x))*(1 - 55*((x)*(x))*(1 - 185*((x)*(x))*(1 - 351*((x)*(x))*(1 - 164*((x)*(x))*(1 - 43*((x)*(x))*(1 - 135*((x)*(x))*(1 - 235*((x)*(x))*(1 - 14*((x)*(x))*(1 - 153*((x)*(x))*(1 - 53*((x)*(x))*(1 - 55*((x)*(x))/351)/150)/253)/15)/171)/68)/15)/39)/55)/18)/3)/5)/3)/16777216;
}

static double legendre_der_29(double _x) {
    long double x = _x;
    return -31556720475*pow(118,(1.0/2.0))*((x)*(x))*(1 - 143*((x)*(x))*(1 - 28*((x)*(x))*(1 - 407*((x)*(x))*(1 - 26*((x)*(x))*(1 - 123*((x)*(x))*(1 - 344*((x)*(x))*(1 - 21*((x)*(x))*(1 - 94*((x)*(x))*(1 - 49*((x)*(x))*(1 - 68*((x)*(x))*(1 - 53*((x)*(x))*(1 - 22*((x)*(x))*(1 - 19*((x)*(x))/126)/65)/92)/77)/38)/51)/8)/91)/22)/3)/28))/2)/33554432 + 145422675*pow(118,(1.0/2.0))/67108864;
}


shape_fn_t legendre_der_tab_1d[] = {

legendre_der_0,
legendre_der_1,
legendre_der_2,
legendre_der_3,
legendre_der_4,
legendre_der_5,
legendre_der_6,
legendre_der_7,
legendre_der_8,
legendre_der_9,
legendre_der_10,
legendre_der_11,
legendre_der_12,
legendre_der_13,
legendre_der_14,
legendre_der_15,
legendre_der_16,
legendre_der_17,
legendre_der_18,
legendre_der_19,
legendre_der_20,
legendre_der_21,
legendre_der_22,
legendre_der_23,
legendre_der_24,
legendre_der_25,
legendre_der_26,
legendre_der_27,
legendre_der_28,
legendre_der_29,
};