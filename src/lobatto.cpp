// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "lobatto.h"

int lobatto_order_1d[] = {
1, 
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


static double lobatto_fn_0(double _x) {
    long double x = _x;
    return 1.0/2.0 - x/2;
}

static double lobatto_fn_1(double _x) {
    long double x = _x;
    return 1.0/2.0 + x/2;
}

static double lobatto_fn_2(double _x) {
    long double x = _x;
    return -pow(6,(1.0/2.0))/4 + pow(6,(1.0/2.0))*((x)*(x))/4;
}

static double lobatto_fn_3(double _x) {
    long double x = _x;
    return -x*pow(10,(1.0/2.0))*(1 - ((x)*(x)))/4;
}

static double lobatto_fn_4(double _x) {
    long double x = _x;
    return -3*pow(14,(1.0/2.0))*((x)*(x))*(1 - 5*((x)*(x))/6)/8 + pow(14,(1.0/2.0))/16;
}

static double lobatto_fn_5(double _x) {
    long double x = _x;
    return 9*x*pow(2,(1.0/2.0))*(1 - 10*((x)*(x))*(1 - 7*((x)*(x))/10)/3)/16;
}

static double lobatto_fn_6(double _x) {
    long double x = _x;
    return 15*pow(22,(1.0/2.0))*((x)*(x))*(1 - 7*((x)*(x))*(1 - 3*((x)*(x))/5)/3)/32 - pow(22,(1.0/2.0))/32;
}

static double lobatto_fn_7(double _x) {
    long double x = _x;
    return -5*x*pow(26,(1.0/2.0))*(1 - 7*((x)*(x))*(1 - 9*((x)*(x))*(1 - 11*((x)*(x))/21)/5))/32;
}

static double lobatto_fn_8(double _x) {
    long double x = _x;
    return -35*pow(30,(1.0/2.0))*((x)*(x))*(1 - 9*((x)*(x))*(1 - 22*((x)*(x))*(1 - 13*((x)*(x))/28)/15)/2)/64 + 5*pow(30,(1.0/2.0))/256;
}

static double lobatto_fn_9(double _x) {
    long double x = _x;
    return 35*x*pow(34,(1.0/2.0))*(1 - 12*((x)*(x))*(1 - 33*((x)*(x))*(1 - 26*((x)*(x))*(1 - 5*((x)*(x))/12)/21)/10))/256;
}

static double lobatto_fn_10(double _x) {
    long double x = _x;
    return 315*pow(38,(1.0/2.0))*((x)*(x))*(1 - 22*((x)*(x))*(1 - 13*((x)*(x))*(1 - 15*((x)*(x))*(1 - 17*((x)*(x))/45)/14)/5)/3)/512 - 7*pow(38,(1.0/2.0))/512;
}

static double lobatto_fn_11(double _x) {
    long double x = _x;
    return -63*x*pow(42,(1.0/2.0))*(1 - 55*((x)*(x))*(1 - 26*((x)*(x))*(1 - 15*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))/55)/18)/7)/5)/3)/512;
}

static double lobatto_fn_12(double _x) {
    long double x = _x;
    return -693*pow(46,(1.0/2.0))*((x)*(x))*(1 - 65*((x)*(x))*(1 - 4*((x)*(x))*(1 - 51*((x)*(x))*(1 - 38*((x)*(x))*(1 - 7*((x)*(x))/22)/45)/28))/6)/1024 + 21*pow(46,(1.0/2.0))/2048;
}

static double lobatto_fn_13(double _x) {
    long double x = _x;
    return 1155*x*pow(2,(1.0/2.0))*(1 - 26*((x)*(x))*(1 - 15*((x)*(x))*(1 - 68*((x)*(x))*(1 - 19*((x)*(x))*(1 - 42*((x)*(x))*(1 - 23*((x)*(x))/78)/55)/12)/21)/2))/2048;
}

static double lobatto_fn_14(double _x) {
    long double x = _x;
    return 9009*pow(6,(1.0/2.0))*((x)*(x))*(1 - 15*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))*(1 - 7*((x)*(x))*(1 - 23*((x)*(x))*(1 - 25*((x)*(x))/91)/33)/5)/7)/3))/4096 - 99*pow(6,(1.0/2.0))/4096;
}

static double lobatto_fn_15(double _x) {
    long double x = _x;
    return -429*x*pow(58,(1.0/2.0))*(1 - 35*((x)*(x))*(1 - 51*((x)*(x))*(1 - 95*((x)*(x))*(1 - 7*((x)*(x))*(1 - 69*((x)*(x))*(1 - 25*((x)*(x))*(1 - 9*((x)*(x))/35)/39)/55)/3)/21)/5))/4096;
}

static double lobatto_fn_16(double _x) {
    long double x = _x;
    return -6435*pow(62,(1.0/2.0))*((x)*(x))*(1 - 119*((x)*(x))*(1 - 38*((x)*(x))*(1 - 15*((x)*(x))*(1 - 92*((x)*(x))*(1 - 25*((x)*(x))*(1 - 54*((x)*(x))*(1 - 29*((x)*(x))/120)/91)/22)/45)/4)/5)/6)/8192 + 429*pow(62,(1.0/2.0))/65536;
}

static double lobatto_fn_17(double _x) {
    long double x = _x;
    return 6435*x*pow(66,(1.0/2.0))*(1 - 136*((x)*(x))*(1 - 133*((x)*(x))*(1 - 6*((x)*(x))*(1 - 115*((x)*(x))*(1 - 20*((x)*(x))*(1 - 27*((x)*(x))*(1 - 58*((x)*(x))*(1 - 31*((x)*(x))/136)/105)/26)/11)/36))/10)/3)/65536;
}

static double lobatto_fn_18(double _x) {
    long double x = _x;
    return 109395*pow(70,(1.0/2.0))*((x)*(x))*(1 - 76*((x)*(x))*(1 - 49*((x)*(x))*(1 - 69*((x)*(x))*(1 - 25*((x)*(x))*(1 - 18*((x)*(x))*(1 - 87*((x)*(x))*(1 - 31*((x)*(x))*(1 - 11*((x)*(x))/51)/60)/91)/11)/9)/14)/5)/3)/131072 - 715*pow(70,(1.0/2.0))/131072;
}

static double lobatto_fn_19(double _x) {
    long double x = _x;
    return -12155*x*pow(74,(1.0/2.0))*(1 - 57*((x)*(x))*(1 - 84*((x)*(x))*(1 - 23*((x)*(x))*(1 - 25*((x)*(x))*(1 - 27*((x)*(x))*(1 - 58*((x)*(x))*(1 - 31*((x)*(x))*(1 - 33*((x)*(x))*(1 - 35*((x)*(x))/171)/68)/35)/39)/11)/6)/3)/5))/131072;
}

static double lobatto_fn_20(double _x) {
    long double x = _x;
    return -230945*pow(78,(1.0/2.0))*((x)*(x))*(1 - 63*((x)*(x))*(1 - 184*((x)*(x))*(1 - 25*((x)*(x))*(1 - 18*((x)*(x))*(1 - 145*((x)*(x))*(1 - 124*((x)*(x))*(1 - 33*((x)*(x))*(1 - 70*((x)*(x))*(1 - 37*((x)*(x))/190)/153)/40)/91)/66)/5)/4)/15)/2)/262144 + 2431*pow(78,(1.0/2.0))/524288;
}

static double lobatto_fn_21(double _x) {
    long double x = _x;
    return 46189*x*pow(82,(1.0/2.0))*(1 - 70*((x)*(x))*(1 - 207*((x)*(x))*(1 - 200*((x)*(x))*(1 - 21*((x)*(x))*(1 - 174*((x)*(x))*(1 - 155*((x)*(x))*(1 - 44*((x)*(x))*(1 - 105*((x)*(x))*(1 - 74*((x)*(x))*(1 - 13*((x)*(x))/70)/171)/136)/35)/78)/55)/4)/21)/10))/524288;
}

static double lobatto_fn_22(double _x) {
    long double x = _x;
    return 969969*pow(86,(1.0/2.0))*((x)*(x))*(1 - 115*((x)*(x))*(1 - 15*((x)*(x))*(1 - 54*((x)*(x))*(1 - 203*((x)*(x))*(1 - 31*((x)*(x))*(1 - 165*((x)*(x))*(1 - 7*((x)*(x))*(1 - 37*((x)*(x))*(1 - 39*((x)*(x))*(1 - 41*((x)*(x))/231)/95)/51)/6)/91)/11)/45)/7))/3)/1048576 - 4199*pow(86,(1.0/2.0))/1048576;
}

static double lobatto_fn_23(double _x) {
    long double x = _x;
    return -264537*x*pow(10,(1.0/2.0))*(1 - 253*((x)*(x))*(1 - 25*((x)*(x))*(1 - 81*((x)*(x))*(1 - 58*((x)*(x))*(1 - 217*((x)*(x))*(1 - 33*((x)*(x))*(1 - 5*((x)*(x))*(1 - 37*((x)*(x))*(1 - 13*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))/253)/105)/19)/34)/3)/13)/55)/9)/7))/3)/1048576;
}

static double lobatto_fn_24(double _x) {
    long double x = _x;
    return -2028117*pow(94,(1.0/2.0))*((x)*(x))*(1 - 275*((x)*(x))*(1 - 18*((x)*(x))*(1 - 261*((x)*(x))*(1 - 248*((x)*(x))*(1 - 7*((x)*(x))*(1 - 30*((x)*(x))*(1 - 37*((x)*(x))*(1 - 52*((x)*(x))*(1 - 123*((x)*(x))*(1 - 86*((x)*(x))*(1 - 15*((x)*(x))/92)/231)/190)/51)/24)/13)/2)/45)/28))/6)/2097152 + 29393*pow(94,(1.0/2.0))/8388608;
}

static double lobatto_fn_25(double _x) {
    long double x = _x;
    return 4732273*x*pow(2,(1.0/2.0))*(1 - 100*((x)*(x))*(1 - 297*((x)*(x))*(1 - 290*((x)*(x))*(1 - 31*((x)*(x))*(1 - 24*((x)*(x))*(1 - 245*((x)*(x))*(1 - 74*((x)*(x))*(1 - 195*((x)*(x))*(1 - 164*((x)*(x))*(1 - 43*((x)*(x))*(1 - 90*((x)*(x))*(1 - 47*((x)*(x))/300)/253)/70)/171)/136)/35)/78)/5)/4)/21)/10))/8388608;
}

static double lobatto_fn_26(double _x) {
    long double x = _x;
    return 16900975*pow(102,(1.0/2.0))*((x)*(x))*(1 - 54*((x)*(x))*(1 - 319*((x)*(x))*(1 - 155*((x)*(x))*(1 - 33*((x)*(x))*(1 - 140*((x)*(x))*(1 - 37*((x)*(x))*(1 - 39*((x)*(x))*(1 - 205*((x)*(x))*(1 - 86*((x)*(x))*(1 - 45*((x)*(x))*(1 - 47*((x)*(x))*(1 - 49*((x)*(x))/325)/138)/77)/95)/153)/20)/13)/33)/5)/14)/15))/16777216 - 52003*pow(102,(1.0/2.0))/16777216;
}

static double lobatto_fn_27(double _x) {
    long double x = _x;
    return -1300075*x*pow(106,(1.0/2.0))*(1 - 117*((x)*(x))*(1 - 174*((x)*(x))*(1 - 341*((x)*(x))*(1 - 55*((x)*(x))*(1 - 63*((x)*(x))*(1 - 148*((x)*(x))*(1 - 13*((x)*(x))*(1 - 123*((x)*(x))*(1 - 215*((x)*(x))*(1 - 6*((x)*(x))*(1 - 141*((x)*(x))*(1 - 49*((x)*(x))*(1 - 17*((x)*(x))/117)/150)/253)/7)/171)/68)/5)/39)/11)/6)/21)/5))/16777216;
}

static double lobatto_fn_28(double _x) {
    long double x = _x;
    return -35102025*pow(110,(1.0/2.0))*((x)*(x))*(1 - 377*((x)*(x))*(1 - 124*((x)*(x))*(1 - 363*((x)*(x))*(1 - 70*((x)*(x))*(1 - 111*((x)*(x))*(1 - 24*((x)*(x))*(1 - 287*((x)*(x))*(1 - 86*((x)*(x))*(1 - 45*((x)*(x))*(1 - 188*((x)*(x))*(1 - 49*((x)*(x))*(1 - 102*((x)*(x))*(1 - 53*((x)*(x))/378)/325)/92)/231)/38)/51)/120)/7)/22)/9)/28)/5)/6)/33554432 + 185725*pow(110,(1.0/2.0))/67108864;
}

static double lobatto_fn_29(double _x) {
    long double x = _x;
    return 5014575*x*pow(114,(1.0/2.0))*(1 - 406*((x)*(x))*(1 - 403*((x)*(x))*(1 - 132*((x)*(x))*(1 - 385*((x)*(x))*(1 - 74*((x)*(x))*(1 - 9*((x)*(x))*(1 - 328*((x)*(x))*(1 - 301*((x)*(x))*(1 - 30*((x)*(x))*(1 - 47*((x)*(x))*(1 - 196*((x)*(x))*(1 - 51*((x)*(x))*(1 - 106*((x)*(x))*(1 - 55*((x)*(x))/406)/351)/100)/253)/42)/19)/136)/105)/2)/11)/36)/7)/10)/3)/67108864;
}


shape_fn_t lobatto_fn_tab_1d[] = {

lobatto_fn_0,
lobatto_fn_1,
lobatto_fn_2,
lobatto_fn_3,
lobatto_fn_4,
lobatto_fn_5,
lobatto_fn_6,
lobatto_fn_7,
lobatto_fn_8,
lobatto_fn_9,
lobatto_fn_10,
lobatto_fn_11,
lobatto_fn_12,
lobatto_fn_13,
lobatto_fn_14,
lobatto_fn_15,
lobatto_fn_16,
lobatto_fn_17,
lobatto_fn_18,
lobatto_fn_19,
lobatto_fn_20,
lobatto_fn_21,
lobatto_fn_22,
lobatto_fn_23,
lobatto_fn_24,
lobatto_fn_25,
lobatto_fn_26,
lobatto_fn_27,
lobatto_fn_28,
lobatto_fn_29,
};



static double lobatto_der_0(double _x) {
    long double x = _x;
    return -1.0/2.0;
}

static double lobatto_der_1(double _x) {
    long double x = _x;
    return 1.0/2.0;
}

static double lobatto_der_2(double _x) {
    long double x = _x;
    return x*pow(6,(1.0/2.0))/2;
}

static double lobatto_der_3(double _x) {
    long double x = _x;
    return -pow(10,(1.0/2.0))/4 + 3*pow(10,(1.0/2.0))*((x)*(x))/4;
}

static double lobatto_der_4(double _x) {
    long double x = _x;
    return -3*x*pow(14,(1.0/2.0))*(1 - 5*((x)*(x))/3)/4;
}

static double lobatto_der_5(double _x) {
    long double x = _x;
    return -45*pow(2,(1.0/2.0))*((x)*(x))*(1 - 7*((x)*(x))/6)/8 + 9*pow(2,(1.0/2.0))/16;
}

static double lobatto_der_6(double _x) {
    long double x = _x;
    return 15*x*pow(22,(1.0/2.0))*(1 - 14*((x)*(x))*(1 - 9*((x)*(x))/10)/3)/16;
}

static double lobatto_der_7(double _x) {
    long double x = _x;
    return 105*pow(26,(1.0/2.0))*((x)*(x))*(1 - 3*((x)*(x))*(1 - 11*((x)*(x))/15))/32 - 5*pow(26,(1.0/2.0))/32;
}

static double lobatto_der_8(double _x) {
    long double x = _x;
    return -35*x*pow(30,(1.0/2.0))*(1 - 9*((x)*(x))*(1 - 11*((x)*(x))*(1 - 13*((x)*(x))/21)/5))/32;
}

static double lobatto_der_9(double _x) {
    long double x = _x;
    return -315*pow(34,(1.0/2.0))*((x)*(x))*(1 - 11*((x)*(x))*(1 - 26*((x)*(x))*(1 - 15*((x)*(x))/28)/15)/2)/64 + 35*pow(34,(1.0/2.0))/256;
}

static double lobatto_der_10(double _x) {
    long double x = _x;
    return 315*x*pow(38,(1.0/2.0))*(1 - 44*((x)*(x))*(1 - 39*((x)*(x))*(1 - 10*((x)*(x))*(1 - 17*((x)*(x))/36)/7)/10)/3)/256;
}

static double lobatto_der_11(double _x) {
    long double x = _x;
    return 3465*pow(42,(1.0/2.0))*((x)*(x))*(1 - 26*((x)*(x))*(1 - 3*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))/45)/14))/3)/512 - 63*pow(42,(1.0/2.0))/512;
}

static double lobatto_der_12(double _x) {
    long double x = _x;
    return -693*x*pow(46,(1.0/2.0))*(1 - 65*((x)*(x))*(1 - 6*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))*(1 - 21*((x)*(x))/55)/18)/7))/3)/512;
}

static double lobatto_der_13(double _x) {
    long double x = _x;
    return -45045*pow(2,(1.0/2.0))*((x)*(x))*(1 - 25*((x)*(x))*(1 - 68*((x)*(x))*(1 - 57*((x)*(x))*(1 - 14*((x)*(x))*(1 - 23*((x)*(x))/66)/15)/28)/15)/2)/1024 + 1155*pow(2,(1.0/2.0))/2048;
}

static double lobatto_der_14(double _x) {
    long double x = _x;
    return 9009*x*pow(6,(1.0/2.0))*(1 - 30*((x)*(x))*(1 - 17*((x)*(x))*(1 - 76*((x)*(x))*(1 - 7*((x)*(x))*(1 - 46*((x)*(x))*(1 - 25*((x)*(x))/78)/55)/4)/21)/2))/2048;
}

static double lobatto_der_15(double _x) {
    long double x = _x;
    return 45045*pow(58,(1.0/2.0))*((x)*(x))*(1 - 17*((x)*(x))*(1 - 19*((x)*(x))*(1 - 3*((x)*(x))*(1 - 23*((x)*(x))*(1 - 25*((x)*(x))*(1 - 27*((x)*(x))/91)/33)/15))/3))/4096 - 429*pow(58,(1.0/2.0))/4096;
}

static double lobatto_der_16(double _x) {
    long double x = _x;
    return -6435*x*pow(62,(1.0/2.0))*(1 - 119*((x)*(x))*(1 - 57*((x)*(x))*(1 - 5*((x)*(x))*(1 - 23*((x)*(x))*(1 - 15*((x)*(x))*(1 - 9*((x)*(x))*(1 - 29*((x)*(x))/105)/13)/11)/9))/5)/3)/4096;
}

static double lobatto_der_17(double _x) {
    long double x = _x;
    return -109395*pow(66,(1.0/2.0))*((x)*(x))*(1 - 133*((x)*(x))*(1 - 42*((x)*(x))*(1 - 115*((x)*(x))*(1 - 20*((x)*(x))*(1 - 27*((x)*(x))*(1 - 58*((x)*(x))*(1 - 31*((x)*(x))/120)/91)/22)/9)/28)/5)/6)/8192 + 6435*pow(66,(1.0/2.0))/65536;
}

static double lobatto_der_18(double _x) {
    long double x = _x;
    return 109395*x*pow(70,(1.0/2.0))*(1 - 152*((x)*(x))*(1 - 147*((x)*(x))*(1 - 46*((x)*(x))*(1 - 125*((x)*(x))*(1 - 108*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 33*((x)*(x))/136)/105)/26)/55)/36)/7)/10)/3)/65536;
}

static double lobatto_der_19(double _x) {
    long double x = _x;
    return 2078505*pow(74,(1.0/2.0))*((x)*(x))*(1 - 28*((x)*(x))*(1 - 161*((x)*(x))*(1 - 75*((x)*(x))*(1 - 3*((x)*(x))*(1 - 58*((x)*(x))*(1 - 93*((x)*(x))*(1 - 11*((x)*(x))*(1 - 35*((x)*(x))/153)/20)/91)/33))/14)/15))/131072 - 12155*pow(74,(1.0/2.0))/131072;
}

static double lobatto_der_20(double _x) {
    long double x = _x;
    return -230945*x*pow(78,(1.0/2.0))*(1 - 63*((x)*(x))*(1 - 92*((x)*(x))*(1 - 25*((x)*(x))*(1 - 9*((x)*(x))*(1 - 29*((x)*(x))*(1 - 62*((x)*(x))*(1 - 33*((x)*(x))*(1 - 35*((x)*(x))*(1 - 37*((x)*(x))/171)/68)/35)/39)/11)/2)/3)/5))/131072;
}

static double lobatto_der_21(double _x) {
    long double x = _x;
    return -4849845*pow(82,(1.0/2.0))*((x)*(x))*(1 - 69*((x)*(x))*(1 - 40*((x)*(x))*(1 - 27*((x)*(x))*(1 - 58*((x)*(x))*(1 - 155*((x)*(x))*(1 - 132*((x)*(x))*(1 - 7*((x)*(x))*(1 - 74*((x)*(x))*(1 - 39*((x)*(x))/190)/153)/8)/91)/66)/15)/4)/3)/2)/262144 + 46189*pow(82,(1.0/2.0))/524288;
}

static double lobatto_der_22(double _x) {
    long double x = _x;
    return 969969*x*pow(86,(1.0/2.0))*(1 - 230*((x)*(x))*(1 - 45*((x)*(x))*(1 - 72*((x)*(x))*(1 - 203*((x)*(x))*(1 - 186*((x)*(x))*(1 - 55*((x)*(x))*(1 - 4*((x)*(x))*(1 - 111*((x)*(x))*(1 - 26*((x)*(x))*(1 - 41*((x)*(x))/210)/57)/136)/3)/26)/55)/36)/7)/2)/3)/524288;
}

static double lobatto_der_23(double _x) {
    long double x = _x;
    return 66927861*pow(10,(1.0/2.0))*((x)*(x))*(1 - 125*((x)*(x))*(1 - 81*((x)*(x))*(1 - 58*((x)*(x))*(1 - 217*((x)*(x))*(1 - 3*((x)*(x))*(1 - 25*((x)*(x))*(1 - 37*((x)*(x))*(1 - 13*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))/231)/95)/17)/30)/13))/45)/7)/5)/3)/1048576 - 264537*pow(10,(1.0/2.0))/1048576;
}

static double lobatto_der_24(double _x) {
    long double x = _x;
    return -2028117*x*pow(94,(1.0/2.0))*(1 - 275*((x)*(x))*(1 - 27*((x)*(x))*(1 - 87*((x)*(x))*(1 - 62*((x)*(x))*(1 - 21*((x)*(x))*(1 - 35*((x)*(x))*(1 - 37*((x)*(x))*(1 - 39*((x)*(x))*(1 - 41*((x)*(x))*(1 - 43*((x)*(x))*(1 - 45*((x)*(x))/253)/105)/57)/34)/21)/13)/5)/9)/7))/3)/1048576;
}

static double lobatto_der_25(double _x) {
    long double x = _x;
    return -354920475*pow(2,(1.0/2.0))*((x)*(x))*(1 - 99*((x)*(x))*(1 - 58*((x)*(x))*(1 - 279*((x)*(x))*(1 - 88*((x)*(x))*(1 - 245*((x)*(x))*(1 - 222*((x)*(x))*(1 - 13*((x)*(x))*(1 - 164*((x)*(x))*(1 - 129*((x)*(x))*(1 - 30*((x)*(x))*(1 - 47*((x)*(x))/276)/77)/190)/153)/8)/91)/66)/15)/28)/3)/2)/2097152 + 4732273*pow(2,(1.0/2.0))/8388608;
}

static double lobatto_der_26(double _x) {
    long double x = _x;
    return 16900975*x*pow(102,(1.0/2.0))*(1 - 108*((x)*(x))*(1 - 319*((x)*(x))*(1 - 310*((x)*(x))*(1 - 33*((x)*(x))*(1 - 56*((x)*(x))*(1 - 259*((x)*(x))*(1 - 78*((x)*(x))*(1 - 205*((x)*(x))*(1 - 172*((x)*(x))*(1 - 9*((x)*(x))*(1 - 94*((x)*(x))*(1 - 49*((x)*(x))/300)/253)/14)/171)/136)/35)/78)/11)/4)/21)/10))/8388608;
}

static double lobatto_der_27(double _x) {
    long double x = _x;
    return 456326325*pow(106,(1.0/2.0))*((x)*(x))*(1 - 58*((x)*(x))*(1 - 341*((x)*(x))*(1 - 165*((x)*(x))*(1 - 7*((x)*(x))*(1 - 148*((x)*(x))*(1 - 3*((x)*(x))*(1 - 41*((x)*(x))*(1 - 215*((x)*(x))*(1 - 18*((x)*(x))*(1 - 47*((x)*(x))*(1 - 49*((x)*(x))*(1 - 51*((x)*(x))/325)/138)/77)/19)/153)/20))/33))/14)/15))/16777216 - 1300075*pow(106,(1.0/2.0))/16777216;
}

static double lobatto_der_28(double _x) {
    long double x = _x;
    return -35102025*x*pow(110,(1.0/2.0))*(1 - 377*((x)*(x))*(1 - 186*((x)*(x))*(1 - 121*((x)*(x))*(1 - 175*((x)*(x))*(1 - 333*((x)*(x))*(1 - 4*((x)*(x))*(1 - 41*((x)*(x))*(1 - 129*((x)*(x))*(1 - 25*((x)*(x))*(1 - 94*((x)*(x))*(1 - 147*((x)*(x))*(1 - 17*((x)*(x))*(1 - 53*((x)*(x))/351)/50)/253)/105)/19)/68)/15))/55)/18)/7)/5)/3)/16777216;
}

static double lobatto_der_29(double _x) {
    long double x = _x;
    return -1017958725*pow(114,(1.0/2.0))*((x)*(x))*(1 - 403*((x)*(x))*(1 - 132*((x)*(x))*(1 - 55*((x)*(x))*(1 - 74*((x)*(x))*(1 - 117*((x)*(x))*(1 - 328*((x)*(x))*(1 - 301*((x)*(x))*(1 - 30*((x)*(x))*(1 - 47*((x)*(x))*(1 - 28*((x)*(x))*(1 - 51*((x)*(x))*(1 - 106*((x)*(x))*(1 - 55*((x)*(x))/378)/325)/92)/33)/38)/17)/120)/91)/22)/9)/4)/5)/6)/33554432 + 5014575*pow(114,(1.0/2.0))/67108864;
}


shape_fn_t lobatto_der_tab_1d[] = {

lobatto_der_0,
lobatto_der_1,
lobatto_der_2,
lobatto_der_3,
lobatto_der_4,
lobatto_der_5,
lobatto_der_6,
lobatto_der_7,
lobatto_der_8,
lobatto_der_9,
lobatto_der_10,
lobatto_der_11,
lobatto_der_12,
lobatto_der_13,
lobatto_der_14,
lobatto_der_15,
lobatto_der_16,
lobatto_der_17,
lobatto_der_18,
lobatto_der_19,
lobatto_der_20,
lobatto_der_21,
lobatto_der_22,
lobatto_der_23,
lobatto_der_24,
lobatto_der_25,
lobatto_der_26,
lobatto_der_27,
lobatto_der_28,
lobatto_der_29,
};