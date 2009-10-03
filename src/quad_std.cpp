// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

// $Id: quad.h 1086 2008-10-21 09:05:44Z jakub $

#ifndef __HERMES2D_QUAD_H
#define __HERMES2D_QUAD_H


/// Quad1D is a base class for all 1D quadrature points.
///
class Quad1D
{
public:
  
  double2* get_points(int order) const { return tables[order]; }
  int get_num_points(int order) const { return np[order]; };
  
  int get_max_order() const { return max_order; }  
  double get_ref_vertex(int n) const { return ref_vert[n]; }

protected:
  
  double2** tables;
  int* np;

  double ref_vert[2];
  int max_order;

  virtual void dummy_fn() = 0; // to prevent this class from being instantiated
};

/// 1D quadrature points on the standard reference domain (-1,1)

class Quad1DStd : public Quad1D
{
  public: Quad1DStd();
    
  virtual void dummy_fn() {}  
};

//// 1D quadrature tables //////////////////////////////

static double2 std_pts_0_1_1d[] = 
{
  { 0.0, 2.0 }
};

static double2 std_pts_2_3_1d[] = 
{
  { -0.57735026918963,  1.0 }, //  { -1.0/sqrt(3.0), 1.0 },
  {  0.57735026918963,  1.0 }  //  {  1.0/sqrt(3.0), 1.0 },
};

static double2 std_pts_4_5_1d[] = 
{
  { -0.77459666924148 /*-sqrt(3.0/5.0)*/,  5.0/9.0 },
  {  0.0,            8.0/9.0 },
  {  0.77459666924148 /*sqrt(3.0/5.0)*/,  5.0/9.0 }
};

static double2 std_pts_6_7_1d[] = 
{
  { -0.86113631159405,  0.34785484513745 },
  { -0.33998104358486,  0.65214515486255 },
  {  0.33998104358486,  0.65214515486255 },
  {  0.86113631159405,  0.34785484513745 }
};

static double2 std_pts_8_9_1d[] = 
{
  { -0.90617984593866,  0.23692688505619 },
  { -0.53846931010568,  0.47862867049937 },
  {  0.00000000000000,  128.0 / 225.0    },
  {  0.53846931010568,  0.47862867049937 },
  {  0.90617984593866,  0.23692688505619 }
};

static double2 std_pts_10_11_1d[] = 
{
  { -0.93246951420315,  0.17132449237917 },
  { -0.66120938646627,  0.36076157304814 },
  { -0.23861918608320,  0.46791393457269 },
  {  0.23861918608320,  0.46791393457269 },
  {  0.66120938646627,  0.36076157304814 },
  {  0.93246951420315,  0.17132449237917 }
};

static double2 std_pts_12_13_1d[] = 
{
  { -0.94910791234276,  0.12948496616887 },
  { -0.74153118559939,  0.27970539148928 },
  { -0.40584515137740,  0.38183005050512 },
  {  0.00000000000000,  0.41795918367347 },
  {  0.40584515137740,  0.38183005050512 },
  {  0.74153118559939,  0.27970539148928 },
  {  0.94910791234276,  0.12948496616887 }
};

static double2 std_pts_14_15_1d[] = 
{
  { -0.96028985649754,  0.10122853629038 },
  { -0.79666647741363,  0.22238103445337 },
  { -0.52553240991633,  0.31370664587789 },
  { -0.18343464249565,  0.36268378337836 },
  {  0.18343464249565,  0.36268378337836 },
  {  0.52553240991633,  0.31370664587789 },
  {  0.79666647741363,  0.22238103445337 },
  {  0.96028985649754,  0.10122853629038 }
};

static double2 std_pts_16_17_1d[] = 
{
  { -0.96816023950763,  0.08127438836157 },
  { -0.83603110732664,  0.18064816069486 },
  { -0.61337143270059,  0.26061069640294 },
  { -0.32425342340381,  0.31234707704000 },
  {  0.00000000000000,  0.33023935500126 },
  {  0.32425342340381,  0.31234707704000 },
  {  0.61337143270059,  0.26061069640294 },
  {  0.83603110732664,  0.18064816069486 },
  {  0.96816023950763,  0.08127438836157 }
};

static double2 std_pts_18_19_1d[] = 
{
  { -0.97390652851717,  0.06667134430869 },
  { -0.86506336668898,  0.14945134915058 },
  { -0.67940956829902,  0.21908636251598 },
  { -0.43339539412925,  0.26926671931000 },
  { -0.14887433898163,  0.29552422471475 },
  {  0.14887433898163,  0.29552422471475 },
  {  0.43339539412925,  0.26926671931000 },
  {  0.67940956829902,  0.21908636251598 },
  {  0.86506336668898,  0.14945134915058 },
  {  0.97390652851717,  0.06667134430869 }
};

static double2 std_pts_20_21_1d[] = 
{
  { -0.97822865814606,  0.05566856711617 },
  { -0.88706259976810,  0.12558036946490 },
  { -0.73015200557405,  0.18629021092773 },
  { -0.51909612920681,  0.23319376459199 },
  { -0.26954315595234,  0.26280454451025 },
  {  0.00000000000000,  0.27292508677790 },
  {  0.26954315595234,  0.26280454451025 },
  {  0.51909612920681,  0.23319376459199 },
  {  0.73015200557405,  0.18629021092773 },
  {  0.88706259976810,  0.12558036946490 },
  {  0.97822865814606,  0.05566856711617 }
};

static double2 std_pts_22_23_1d[] = 
{
  { -0.98156063424672,  0.04717533638651 },
  { -0.90411725637047,  0.10693932599532 },
  { -0.76990267419430,  0.16007832854335 },
  { -0.58731795428662,  0.20316742672307 },
  { -0.36783149899818,  0.23349253653835 },
  { -0.12523340851147,  0.24914704581340 },
  {  0.12523340851147,  0.24914704581340 },
  {  0.36783149899818,  0.23349253653835 },
  {  0.58731795428662,  0.20316742672307 },
  {  0.76990267419430,  0.16007832854335 },
  {  0.90411725637047,  0.10693932599532 },
  {  0.98156063424672,  0.04717533638651 }
};

static double2 std_pts_24_1d[] = 
{
  { -0.98418305471859,  0.04048400476532 },
  { -0.91759839922298,  0.09212149983773 },
  { -0.80157809073331,  0.13887351021979 },
  { -0.64234933944034,  0.17814598076195 },
  { -0.44849275103645,  0.20781604753689 },
  { -0.23045831595513,  0.22628318026290 },
  {  0.00000000000000,  0.23255155323087 },
  {  0.23045831595513,  0.22628318026290 },
  {  0.44849275103645,  0.20781604753689 },
  {  0.64234933944034,  0.17814598076195 },
  {  0.80157809073331,  0.13887351021979 },
  {  0.91759839922298,  0.09212149983773 },
  {  0.98418305471859,  0.04048400476532 }
};



static double2* std_tables_1d[] =
{
  std_pts_0_1_1d,   std_pts_0_1_1d,
  std_pts_2_3_1d,   std_pts_2_3_1d,
  std_pts_4_5_1d,   std_pts_4_5_1d,
  std_pts_6_7_1d,   std_pts_6_7_1d,
  std_pts_8_9_1d,   std_pts_8_9_1d,
  std_pts_10_11_1d, std_pts_10_11_1d,
  std_pts_12_13_1d, std_pts_12_13_1d,
  std_pts_14_15_1d, std_pts_14_15_1d,
  std_pts_16_17_1d, std_pts_16_17_1d,
  std_pts_18_19_1d, std_pts_18_19_1d,
  std_pts_20_21_1d, std_pts_20_21_1d,
  std_pts_22_23_1d, std_pts_22_23_1d,
  std_pts_24_1d
};


static int std_np_1d[] =
{
  sizeof(std_pts_0_1_1d) / sizeof(double2),
  sizeof(std_pts_0_1_1d) / sizeof(double2),
  sizeof(std_pts_2_3_1d) / sizeof(double2),
  sizeof(std_pts_2_3_1d) / sizeof(double2),
  sizeof(std_pts_4_5_1d) / sizeof(double2),
  sizeof(std_pts_4_5_1d) / sizeof(double2),
  sizeof(std_pts_6_7_1d) / sizeof(double2),
  sizeof(std_pts_6_7_1d) / sizeof(double2),
  sizeof(std_pts_8_9_1d) / sizeof(double2),
  sizeof(std_pts_8_9_1d) / sizeof(double2),
  sizeof(std_pts_10_11_1d) / sizeof(double2),
  sizeof(std_pts_10_11_1d) / sizeof(double2),
  sizeof(std_pts_12_13_1d) / sizeof(double2),
  sizeof(std_pts_12_13_1d) / sizeof(double2),
  sizeof(std_pts_14_15_1d) / sizeof(double2),
  sizeof(std_pts_14_15_1d) / sizeof(double2),
  sizeof(std_pts_16_17_1d) / sizeof(double2),
  sizeof(std_pts_16_17_1d) / sizeof(double2),
  sizeof(std_pts_18_19_1d) / sizeof(double2),
  sizeof(std_pts_18_19_1d) / sizeof(double2),
  sizeof(std_pts_20_21_1d) / sizeof(double2),
  sizeof(std_pts_20_21_1d) / sizeof(double2),
  sizeof(std_pts_22_23_1d) / sizeof(double2),
  sizeof(std_pts_22_23_1d) / sizeof(double2),
  sizeof(std_pts_24_1d) / sizeof(double2)
};


Quad1DStd::Quad1DStd()
{
  tables = std_tables_1d;
  np = std_np_1d;
  ref_vert[0] = -1.0;
  ref_vert[1] = 1.0;
  max_order = 24;
}

Quad1DStd g_quad_1d_std;

#endif
