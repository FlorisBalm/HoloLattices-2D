num_EOMs=9

colPositions_NoDOF =\
"""
                    col[0].k = k + -2; col[0].j = j + -2; col[0].i = i + 0;
                    col[1].k = k + -2; col[1].j = j + -1; col[1].i = i + 0;
                    col[2].k = k + -2; col[2].j = j + 0; col[2].i = i + -2;
                    col[3].k = k + -2; col[3].j = j + 0; col[3].i = i + -1;
                    col[4].k = k + -2; col[4].j = j + 0; col[4].i = i + 0;
                    col[5].k = k + -2; col[5].j = j + 0; col[5].i = i + 1;
                    col[6].k = k + -2; col[6].j = j + 0; col[6].i = i + 2;
                    col[7].k = k + -2; col[7].j = j + 1; col[7].i = i + 0;
                    col[8].k = k + -2; col[8].j = j + 2; col[8].i = i + 0;
                    col[9].k = k + -1; col[9].j = j + -2; col[9].i = i + 0;
                    col[10].k = k + -1; col[10].j = j + -1; col[10].i = i + 0;
                    col[11].k = k + -1; col[11].j = j + 0; col[11].i = i + -2;
                    col[12].k = k + -1; col[12].j = j + 0; col[12].i = i + -1;
                    col[13].k = k + -1; col[13].j = j + 0; col[13].i = i + 0;
                    col[14].k = k + -1; col[14].j = j + 0; col[14].i = i + 1;
                    col[15].k = k + -1; col[15].j = j + 0; col[15].i = i + 2;
                    col[16].k = k + -1; col[16].j = j + 1; col[16].i = i + 0;
                    col[17].k = k + -1; col[17].j = j + 2; col[17].i = i + 0;
                    col[18].k = k + 0; col[18].j = j + -2; col[18].i = i + -2;
                    col[19].k = k + 0; col[19].j = j + -2; col[19].i = i + -1;
                    col[20].k = k + 0; col[20].j = j + -2; col[20].i = i + 0;
                    col[21].k = k + 0; col[21].j = j + -2; col[21].i = i + 1;
                    col[22].k = k + 0; col[22].j = j + -2; col[22].i = i + 2;
                    col[23].k = k + 0; col[23].j = j + -1; col[23].i = i + -2;
                    col[24].k = k + 0; col[24].j = j + -1; col[24].i = i + -1;
                    col[25].k = k + 0; col[25].j = j + -1; col[25].i = i + 0;
                    col[26].k = k + 0; col[26].j = j + -1; col[26].i = i + 1;
                    col[27].k = k + 0; col[27].j = j + -1; col[27].i = i + 2;
                    col[28].k = k + 0; col[28].j = j + 0; col[28].i = i + -2;
                    col[29].k = k + 0; col[29].j = j + 0; col[29].i = i + -1;
                    col[30].k = k + 0; col[30].j = j + 0; col[30].i = i + 0;
                    col[31].k = k + 0; col[31].j = j + 0; col[31].i = i + 1;
                    col[32].k = k + 0; col[32].j = j + 0; col[32].i = i + 2;
                    col[33].k = k + 0; col[33].j = j + 1; col[33].i = i + -2;
                    col[34].k = k + 0; col[34].j = j + 1; col[34].i = i + -1;
                    col[35].k = k + 0; col[35].j = j + 1; col[35].i = i + 0;
                    col[36].k = k + 0; col[36].j = j + 1; col[36].i = i + 1;
                    col[37].k = k + 0; col[37].j = j + 1; col[37].i = i + 2;
                    col[38].k = k + 0; col[38].j = j + 2; col[38].i = i + -2;
                    col[39].k = k + 0; col[39].j = j + 2; col[39].i = i + -1;
                    col[40].k = k + 0; col[40].j = j + 2; col[40].i = i + 0;
                    col[41].k = k + 0; col[41].j = j + 2; col[41].i = i + 1;
                    col[42].k = k + 0; col[42].j = j + 2; col[42].i = i + 2;
                    col[43].k = k + 1; col[43].j = j + -2; col[43].i = i + 0;
                    col[44].k = k + 1; col[44].j = j + -1; col[44].i = i + 0;
                    col[45].k = k + 1; col[45].j = j + 0; col[45].i = i + -2;
                    col[46].k = k + 1; col[46].j = j + 0; col[46].i = i + -1;
                    col[47].k = k + 1; col[47].j = j + 0; col[47].i = i + 0;
                    col[48].k = k + 1; col[48].j = j + 0; col[48].i = i + 1;
                    col[49].k = k + 1; col[49].j = j + 0; col[49].i = i + 2;
                    col[50].k = k + 1; col[50].j = j + 1; col[50].i = i + 0;
                    col[51].k = k + 1; col[51].j = j + 2; col[51].i = i + 0;
                    col[52].k = k + 2; col[52].j = j + -2; col[52].i = i + 0;
                    col[53].k = k + 2; col[53].j = j + -1; col[53].i = i + 0;
                    col[54].k = k + 2; col[54].j = j + 0; col[54].i = i + -2;
                    col[55].k = k + 2; col[55].j = j + 0; col[55].i = i + -1;
                    col[56].k = k + 2; col[56].j = j + 0; col[56].i = i + 0;
                    col[57].k = k + 2; col[57].j = j + 0; col[57].i = i + 1;
                    col[58].k = k + 2; col[58].j = j + 0; col[58].i = i + 2;
                    col[59].k = k + 2; col[59].j = j + 1; col[59].i = i + 0;
                    col[60].k = k + 2; col[60].j = j + 2; col[60].i = i + 0;

"""


colPositions_DOFOnly =\
        """
for(int kkk = 0; kkk < 61; ++kkk){{
    col[kkk].c = {n};
}}
"""

valExpression_bulk = """

                    row.c = {eom};

                    cVal[0] = {C}0{E};
                    cVal[1] = {C}1{E};
                    cVal[2] = {C}2{E};
                    cVal[3] = {C}3{E};
                    cVal[4] = {C}4{E};
                    cVal[5] = {C}5{E};
                    cVal[6] = {C}6{E};
                    cVal[7] = {C}7{E};
                    cVal[8] = {C}8{E};
                    cVal[9] = {C}9{E};

                    val[0] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[1] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[2] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[3] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[4] = 1.*cVal[3]*derivativeCoefsZ[1][0] + 1.*cVal[6]*derivativeCoefsZ[2][0];
                    val[5] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[6] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[7] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[8] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[9] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[10] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[11] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[12] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[13] = 1.*cVal[3]*derivativeCoefsZ[1][1] + 1.*cVal[6]*derivativeCoefsZ[2][1];
                    val[14] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[15] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[16] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[17] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[18] = (0.006944444444444444*cVal[7])/(dx*dy);
                    val[19] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[20] = (0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) + (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[21] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[22] = (-0.006944444444444444*cVal[7])/(dx*dy);
                    val[23] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[24] = (0.4444444444444444*cVal[7])/(dx*dy);
                    val[25] = (-0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) - (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[26] = (-0.4444444444444444*cVal[7])/(dx*dy);
                    val[27] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[28] = (0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) + (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[29] = (-0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) - (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[30] = 1.*cVal[0] - (2.5*cVal[4])/(dx*dx) - (2.5*cVal[5])/(dy*dy) + 1.*cVal[3]*derivativeCoefsZ[1][2] + 1.*cVal[6]*derivativeCoefsZ[2][2];
                    val[31] = (0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) + (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[32] = (-0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) - (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[33] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[34] = (-0.4444444444444444*cVal[7])/(dx*dy);
                    val[35] = (0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) + (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[36] = (0.4444444444444444*cVal[7])/(dx*dy);
                    val[37] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[38] = (-0.006944444444444444*cVal[7])/(dx*dy);
                    val[39] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[40] = (-0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) - (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[41] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[42] = (0.006944444444444444*cVal[7])/(dx*dy);
                    val[43] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[44] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[45] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[46] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[47] = 1.*cVal[3]*derivativeCoefsZ[1][3] + 1.*cVal[6]*derivativeCoefsZ[2][3];
                    val[48] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[49] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[50] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[51] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[52] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][4])/dy;
                    val[53] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][4])/dy;
                    val[54] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][4])/dx;
                    val[55] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][4])/dx;
                    val[56] = 1.*cVal[3]*derivativeCoefsZ[1][4] + 1.*cVal[6]*derivativeCoefsZ[2][4];
                    val[57] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][4])/dx;
                    val[58] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][4])/dx;
                    val[59] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][4])/dy;
                    val[60] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][4])/dy;
                    MatSetValuesStencil(jac, 1, &row, 61, col, val, INSERT_VALUES); 
        """.format(eom="{eom}",C="Coefs_I_{eom}_", E="_{field}(derivs, mu,Q,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")

f = open("HeadersSecondOrder/ExpressionsCoefficientsInternal.h", "w")
#f.write("""
#        #ifndef __EXP_COEF_INT_H
#        #define __EXP_COEF_INT_H
#        """)
f.write(colPositions_NoDOF)
f.write("\n")

for n_field in range(num_EOMs):
    f.write(colPositions_DOFOnly.format(n=n_field))

    for n_eom in range(num_EOMs):
        f.write(valExpression_bulk.format(field=n_field, eom=n_eom))

f.write("\n")
f.close()
##########################################
####
####
####
#### ONE FROM BDY at ZERO SIDE
####
####
####
##########################################

colPositions_NoDOF =\
"""
                    col[0].k = k + -1; col[0].j = j + -2; col[0].i = i + 0;
                    col[1].k = k + -1; col[1].j = j + -1; col[1].i = i + 0;
                    col[2].k = k + -1; col[2].j = j + 0; col[2].i = i + -2;
                    col[3].k = k + -1; col[3].j = j + 0; col[3].i = i + -1;
                    col[4].k = k + -1; col[4].j = j + 0; col[4].i = i + 0;
                    col[5].k = k + -1; col[5].j = j + 0; col[5].i = i + 1;
                    col[6].k = k + -1; col[6].j = j + 0; col[6].i = i + 2;
                    col[7].k = k + -1; col[7].j = j + 1; col[7].i = i + 0;
                    col[8].k = k + -1; col[8].j = j + 2; col[8].i = i + 0;
                    col[9].k = k + 0; col[9].j = j + -2; col[9].i = i + -2;
                    col[10].k = k + 0; col[10].j = j + -2; col[10].i = i + -1;
                    col[11].k = k + 0; col[11].j = j + -2; col[11].i = i + 0;
                    col[12].k = k + 0; col[12].j = j + -2; col[12].i = i + 1;
                    col[13].k = k + 0; col[13].j = j + -2; col[13].i = i + 2;
                    col[14].k = k + 0; col[14].j = j + -1; col[14].i = i + -2;
                    col[15].k = k + 0; col[15].j = j + -1; col[15].i = i + -1;
                    col[16].k = k + 0; col[16].j = j + -1; col[16].i = i + 0;
                    col[17].k = k + 0; col[17].j = j + -1; col[17].i = i + 1;
                    col[18].k = k + 0; col[18].j = j + -1; col[18].i = i + 2;
                    col[19].k = k + 0; col[19].j = j + 0; col[19].i = i + -2;
                    col[20].k = k + 0; col[20].j = j + 0; col[20].i = i + -1;
                    col[21].k = k + 0; col[21].j = j + 0; col[21].i = i + 0;
                    col[22].k = k + 0; col[22].j = j + 0; col[22].i = i + 1;
                    col[23].k = k + 0; col[23].j = j + 0; col[23].i = i + 2;
                    col[24].k = k + 0; col[24].j = j + 1; col[24].i = i + -2;
                    col[25].k = k + 0; col[25].j = j + 1; col[25].i = i + -1;
                    col[26].k = k + 0; col[26].j = j + 1; col[26].i = i + 0;
                    col[27].k = k + 0; col[27].j = j + 1; col[27].i = i + 1;
                    col[28].k = k + 0; col[28].j = j + 1; col[28].i = i + 2;
                    col[29].k = k + 0; col[29].j = j + 2; col[29].i = i + -2;
                    col[30].k = k + 0; col[30].j = j + 2; col[30].i = i + -1;
                    col[31].k = k + 0; col[31].j = j + 2; col[31].i = i + 0;
                    col[32].k = k + 0; col[32].j = j + 2; col[32].i = i + 1;
                    col[33].k = k + 0; col[33].j = j + 2; col[33].i = i + 2;
                    col[34].k = k + 1; col[34].j = j + -2; col[34].i = i + 0;
                    col[35].k = k + 1; col[35].j = j + -1; col[35].i = i + 0;
                    col[36].k = k + 1; col[36].j = j + 0; col[36].i = i + -2;
                    col[37].k = k + 1; col[37].j = j + 0; col[37].i = i + -1;
                    col[38].k = k + 1; col[38].j = j + 0; col[38].i = i + 0;
                    col[39].k = k + 1; col[39].j = j + 0; col[39].i = i + 1;
                    col[40].k = k + 1; col[40].j = j + 0; col[40].i = i + 2;
                    col[41].k = k + 1; col[41].j = j + 1; col[41].i = i + 0;
                    col[42].k = k + 1; col[42].j = j + 2; col[42].i = i + 0;
                    col[43].k = k + 2; col[43].j = j + -2; col[43].i = i + 0;
                    col[44].k = k + 2; col[44].j = j + -1; col[44].i = i + 0;
                    col[45].k = k + 2; col[45].j = j + 0; col[45].i = i + -2;
                    col[46].k = k + 2; col[46].j = j + 0; col[46].i = i + -1;
                    col[47].k = k + 2; col[47].j = j + 0; col[47].i = i + 0;
                    col[48].k = k + 2; col[48].j = j + 0; col[48].i = i + 1;
                    col[49].k = k + 2; col[49].j = j + 0; col[49].i = i + 2;
                    col[50].k = k + 2; col[50].j = j + 1; col[50].i = i + 0;
                    col[51].k = k + 2; col[51].j = j + 2; col[51].i = i + 0;

"""

colPositions_DOFOnly =\
        """
for(int kkk = 0; kkk < 52; ++kkk){{
    col[kkk].c = {n};
}}
"""

valExpression_bulk = """

                    row.c = {eom};

                    cVal[0] = {C}0{E};
                    cVal[1] = {C}1{E};
                    cVal[2] = {C}2{E};
                    cVal[3] = {C}3{E};
                    cVal[4] = {C}4{E};
                    cVal[5] = {C}5{E};
                    cVal[6] = {C}6{E};
                    cVal[7] = {C}7{E};
                    cVal[8] = {C}8{E};
                    cVal[9] = {C}9{E};

                    val[0] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[1] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[2] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[3] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[4] = 1.*cVal[3]*derivativeCoefsZ[1][0] + 1.*cVal[6]*derivativeCoefsZ[2][0];
                    val[5] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[6] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[7] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[8] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[9] = (0.006944444444444444*cVal[7])/(dx*dy);
                    val[10] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[11] = (0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) + (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[12] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[13] = (-0.006944444444444444*cVal[7])/(dx*dy);
                    val[14] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[15] = (0.4444444444444444*cVal[7])/(dx*dy);
                    val[16] = (-0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) - (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[17] = (-0.4444444444444444*cVal[7])/(dx*dy);
                    val[18] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[19] = (0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) + (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[20] = (-0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) - (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[21] = 1.*cVal[0] - (2.5*cVal[4])/(dx*dx) - (2.5*cVal[5])/(dy*dy) + 1.*cVal[3]*derivativeCoefsZ[1][1] + 1.*cVal[6]*derivativeCoefsZ[2][1];
                    val[22] = (0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) + (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[23] = (-0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) - (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[24] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[25] = (-0.4444444444444444*cVal[7])/(dx*dy);
                    val[26] = (0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) + (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[27] = (0.4444444444444444*cVal[7])/(dx*dy);
                    val[28] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[29] = (-0.006944444444444444*cVal[7])/(dx*dy);
                    val[30] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[31] = (-0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) - (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[32] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[33] = (0.006944444444444444*cVal[7])/(dx*dy);
                    val[34] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[35] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[36] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[37] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[38] = 1.*cVal[3]*derivativeCoefsZ[1][2] + 1.*cVal[6]*derivativeCoefsZ[2][2];
                    val[39] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[40] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[41] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[42] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[43] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[44] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[45] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[46] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[47] = 1.*cVal[3]*derivativeCoefsZ[1][3] + 1.*cVal[6]*derivativeCoefsZ[2][3];
                    val[48] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[49] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[50] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[51] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][3])/dy;

                    MatSetValuesStencil(jac, 1, &row, 52, col, val, INSERT_VALUES); 
        """.format(eom="{eom}",C="Coefs_I_{eom}_", E="_{field}(derivs, mu,Q,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")

f = open("HeadersSecondOrder/ExpressionCoefficients_2_0_1.h", "w")
#f.write("""
#        #ifndef __EXP_COEF_INT_H
#        #define __EXP_COEF_INT_H
#        """)
f.write(colPositions_NoDOF)
f.write("\n")

for n_field in range(num_EOMs):
    f.write(colPositions_DOFOnly.format(n=n_field))

    for n_eom in range(num_EOMs):
        f.write(valExpression_bulk.format(field=n_field, eom=n_eom))

f.write("\n")
f.close()


##########################################
####
####
####
#### ONE FROM BDY at ONE SIDE
####
####
####
##########################################

colPositions_NoDOF =\
"""
                    col[0].k = k + -2; col[0].j = j + -2; col[0].i = i + 0;
                    col[1].k = k + -2; col[1].j = j + -1; col[1].i = i + 0;
                    col[2].k = k + -2; col[2].j = j + 0; col[2].i = i + -2;
                    col[3].k = k + -2; col[3].j = j + 0; col[3].i = i + -1;
                    col[4].k = k + -2; col[4].j = j + 0; col[4].i = i + 0;
                    col[5].k = k + -2; col[5].j = j + 0; col[5].i = i + 1;
                    col[6].k = k + -2; col[6].j = j + 0; col[6].i = i + 2;
                    col[7].k = k + -2; col[7].j = j + 1; col[7].i = i + 0;
                    col[8].k = k + -2; col[8].j = j + 2; col[8].i = i + 0;
                    col[9].k = k + -1; col[9].j = j + -2; col[9].i = i + 0;
                    col[10].k = k + -1; col[10].j = j + -1; col[10].i = i + 0;
                    col[11].k = k + -1; col[11].j = j + 0; col[11].i = i + -2;
                    col[12].k = k + -1; col[12].j = j + 0; col[12].i = i + -1;
                    col[13].k = k + -1; col[13].j = j + 0; col[13].i = i + 0;
                    col[14].k = k + -1; col[14].j = j + 0; col[14].i = i + 1;
                    col[15].k = k + -1; col[15].j = j + 0; col[15].i = i + 2;
                    col[16].k = k + -1; col[16].j = j + 1; col[16].i = i + 0;
                    col[17].k = k + -1; col[17].j = j + 2; col[17].i = i + 0;
                    col[18].k = k + 0; col[18].j = j + -2; col[18].i = i + -2;
                    col[19].k = k + 0; col[19].j = j + -2; col[19].i = i + -1;
                    col[20].k = k + 0; col[20].j = j + -2; col[20].i = i + 0;
                    col[21].k = k + 0; col[21].j = j + -2; col[21].i = i + 1;
                    col[22].k = k + 0; col[22].j = j + -2; col[22].i = i + 2;
                    col[23].k = k + 0; col[23].j = j + -1; col[23].i = i + -2;
                    col[24].k = k + 0; col[24].j = j + -1; col[24].i = i + -1;
                    col[25].k = k + 0; col[25].j = j + -1; col[25].i = i + 0;
                    col[26].k = k + 0; col[26].j = j + -1; col[26].i = i + 1;
                    col[27].k = k + 0; col[27].j = j + -1; col[27].i = i + 2;
                    col[28].k = k + 0; col[28].j = j + 0; col[28].i = i + -2;
                    col[29].k = k + 0; col[29].j = j + 0; col[29].i = i + -1;
                    col[30].k = k + 0; col[30].j = j + 0; col[30].i = i + 0;
                    col[31].k = k + 0; col[31].j = j + 0; col[31].i = i + 1;
                    col[32].k = k + 0; col[32].j = j + 0; col[32].i = i + 2;
                    col[33].k = k + 0; col[33].j = j + 1; col[33].i = i + -2;
                    col[34].k = k + 0; col[34].j = j + 1; col[34].i = i + -1;
                    col[35].k = k + 0; col[35].j = j + 1; col[35].i = i + 0;
                    col[36].k = k + 0; col[36].j = j + 1; col[36].i = i + 1;
                    col[37].k = k + 0; col[37].j = j + 1; col[37].i = i + 2;
                    col[38].k = k + 0; col[38].j = j + 2; col[38].i = i + -2;
                    col[39].k = k + 0; col[39].j = j + 2; col[39].i = i + -1;
                    col[40].k = k + 0; col[40].j = j + 2; col[40].i = i + 0;
                    col[41].k = k + 0; col[41].j = j + 2; col[41].i = i + 1;
                    col[42].k = k + 0; col[42].j = j + 2; col[42].i = i + 2;
                    col[43].k = k + 1; col[43].j = j + -2; col[43].i = i + 0;
                    col[44].k = k + 1; col[44].j = j + -1; col[44].i = i + 0;
                    col[45].k = k + 1; col[45].j = j + 0; col[45].i = i + -2;
                    col[46].k = k + 1; col[46].j = j + 0; col[46].i = i + -1;
                    col[47].k = k + 1; col[47].j = j + 0; col[47].i = i + 0;
                    col[48].k = k + 1; col[48].j = j + 0; col[48].i = i + 1;
                    col[49].k = k + 1; col[49].j = j + 0; col[49].i = i + 2;
                    col[50].k = k + 1; col[50].j = j + 1; col[50].i = i + 0;
                    col[51].k = k + 1; col[51].j = j + 2; col[51].i = i + 0;    

"""

colPositions_DOFOnly =\
        """
for(int kkk = 0; kkk < 52; ++kkk){{
    col[kkk].c = {n};
}}
"""

valExpression_bulk = """

                    row.c = {eom};

                    cVal[0] = {C}0{E};
                    cVal[1] = {C}1{E};
                    cVal[2] = {C}2{E};
                    cVal[3] = {C}3{E};
                    cVal[4] = {C}4{E};
                    cVal[5] = {C}5{E};
                    cVal[6] = {C}6{E};
                    cVal[7] = {C}7{E};
                    cVal[8] = {C}8{E};
                    cVal[9] = {C}9{E};

                    val[0] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[1] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[2] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[3] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[4] = 1.*cVal[3]*derivativeCoefsZ[1][0] + 1.*cVal[6]*derivativeCoefsZ[2][0];
                    val[5] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[6] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
                    val[7] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[8] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
                    val[9] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[10] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[11] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[12] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[13] = 1.*cVal[3]*derivativeCoefsZ[1][1] + 1.*cVal[6]*derivativeCoefsZ[2][1];
                    val[14] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[15] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
                    val[16] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[17] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
                    val[18] = (0.006944444444444444*cVal[7])/(dx*dy);
                    val[19] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[20] = (0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) + (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[21] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[22] = (-0.006944444444444444*cVal[7])/(dx*dy);
                    val[23] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[24] = (0.4444444444444444*cVal[7])/(dx*dy);
                    val[25] = (-0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) - (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[26] = (-0.4444444444444444*cVal[7])/(dx*dy);
                    val[27] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[28] = (0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) + (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[29] = (-0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) - (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[30] = 1.*cVal[0] - (2.5*cVal[4])/(dx*dx) - (2.5*cVal[5])/(dy*dy) + 1.*cVal[3]*derivativeCoefsZ[1][2] + 1.*cVal[6]*derivativeCoefsZ[2][2];
                    val[31] = (0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) + (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[32] = (-0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) - (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
                    val[33] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[34] = (-0.4444444444444444*cVal[7])/(dx*dy);
                    val[35] = (0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) + (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[36] = (0.4444444444444444*cVal[7])/(dx*dy);
                    val[37] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[38] = (-0.006944444444444444*cVal[7])/(dx*dy);
                    val[39] = (0.05555555555555555*cVal[7])/(dx*dy);
                    val[40] = (-0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) - (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
                    val[41] = (-0.05555555555555555*cVal[7])/(dx*dy);
                    val[42] = (0.006944444444444444*cVal[7])/(dx*dy);
                    val[43] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[44] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[45] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[46] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[47] = 1.*cVal[3]*derivativeCoefsZ[1][3] + 1.*cVal[6]*derivativeCoefsZ[2][3];
                    val[48] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[49] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][3])/dx;
                    val[50] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][3])/dy;
                    val[51] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][3])/dy;

                    MatSetValuesStencil(jac, 1, &row, 52, col, val, INSERT_VALUES); 
        """.format(eom="{eom}",C="Coefs_I_{eom}_", E="_{field}(derivs, mu,Q,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")

f = open("HeadersSecondOrder/ExpressionCoefficients_2_1_1.h", "w")
#f.write("""
#        #ifndef __EXP_COEF_INT_H
#        #define __EXP_COEF_INT_H
#        """)
f.write(colPositions_NoDOF)
f.write("\n")

for n_field in range(num_EOMs):
    f.write(colPositions_DOFOnly.format(n=n_field))

    for n_eom in range(num_EOMs):
        f.write(valExpression_bulk.format(field=n_field, eom=n_eom))

f.write("\n")
f.close()
########################
####
####
#### BOUNDARIES
####
####
########################

formatExpression_2_0_NoDOF =\
"""
                    col[0].k = k + 0; col[0].j = j + -2; col[0].i = i + -2;
                    col[1].k = k + 0; col[1].j = j + -2; col[1].i = i + -1;
                    col[2].k = k + 0; col[2].j = j + -2; col[2].i = i + 0;
                    col[3].k = k + 0; col[3].j = j + -2; col[3].i = i + 1;
                    col[4].k = k + 0; col[4].j = j + -2; col[4].i = i + 2;
                    col[5].k = k + 0; col[5].j = j + -1; col[5].i = i + -2;
                    col[6].k = k + 0; col[6].j = j + -1; col[6].i = i + -1;
                    col[7].k = k + 0; col[7].j = j + -1; col[7].i = i + 0;
                    col[8].k = k + 0; col[8].j = j + -1; col[8].i = i + 1;
                    col[9].k = k + 0; col[9].j = j + -1; col[9].i = i + 2;
                    col[10].k = k + 0; col[10].j = j + 0; col[10].i = i + -2;
                    col[11].k = k + 0; col[11].j = j + 0; col[11].i = i + -1;
                    col[12].k = k + 0; col[12].j = j + 0; col[12].i = i + 0;
                    col[13].k = k + 0; col[13].j = j + 0; col[13].i = i + 1;
                    col[14].k = k + 0; col[14].j = j + 0; col[14].i = i + 2;
                    col[15].k = k + 0; col[15].j = j + 1; col[15].i = i + -2;
                    col[16].k = k + 0; col[16].j = j + 1; col[16].i = i + -1;
                    col[17].k = k + 0; col[17].j = j + 1; col[17].i = i + 0;
                    col[18].k = k + 0; col[18].j = j + 1; col[18].i = i + 1;
                    col[19].k = k + 0; col[19].j = j + 1; col[19].i = i + 2;
                    col[20].k = k + 0; col[20].j = j + 2; col[20].i = i + -2;
                    col[21].k = k + 0; col[21].j = j + 2; col[21].i = i + -1;
                    col[22].k = k + 0; col[22].j = j + 2; col[22].i = i + 0;
                    col[23].k = k + 0; col[23].j = j + 2; col[23].i = i + 1;
                    col[24].k = k + 0; col[24].j = j + 2; col[24].i = i + 2;
                    col[25].k = k + 1; col[25].j = j + -2; col[25].i = i + 0;
                    col[26].k = k + 1; col[26].j = j + -1; col[26].i = i + 0;
                    col[27].k = k + 1; col[27].j = j + 0; col[27].i = i + -2;
                    col[28].k = k + 1; col[28].j = j + 0; col[28].i = i + -1;
                    col[29].k = k + 1; col[29].j = j + 0; col[29].i = i + 0;
                    col[30].k = k + 1; col[30].j = j + 0; col[30].i = i + 1;
                    col[31].k = k + 1; col[31].j = j + 0; col[31].i = i + 2;
                    col[32].k = k + 1; col[32].j = j + 1; col[32].i = i + 0;
                    col[33].k = k + 1; col[33].j = j + 2; col[33].i = i + 0;
                    col[34].k = k + 2; col[34].j = j + -2; col[34].i = i + 0;
                    col[35].k = k + 2; col[35].j = j + -1; col[35].i = i + 0;
                    col[36].k = k + 2; col[36].j = j + 0; col[36].i = i + -2;
                    col[37].k = k + 2; col[37].j = j + 0; col[37].i = i + -1;
                    col[38].k = k + 2; col[38].j = j + 0; col[38].i = i + 0;
                    col[39].k = k + 2; col[39].j = j + 0; col[39].i = i + 1;
                    col[40].k = k + 2; col[40].j = j + 0; col[40].i = i + 2;
                    col[41].k = k + 2; col[41].j = j + 1; col[41].i = i + 0;
                    col[42].k = k + 2; col[42].j = j + 2; col[42].i = i + 0;

"""

formatExpression_Bdy_DOFOnly = \
        """
for(int kkk = 0; kkk < 43; ++kkk){{
    col[kkk].c = {n};
}}
"""
valExpression_2_0=\
        """
        row.c = {eom};

        cVal[0] = {C}0{E};
        cVal[1] = {C}1{E};
        cVal[2] = {C}2{E};
        cVal[3] = {C}3{E};
        cVal[4] = {C}4{E};
        cVal[5] = {C}5{E};
        cVal[6] = {C}6{E};
        cVal[7] = {C}7{E};
        cVal[8] = {C}8{E};
        cVal[9] = {C}9{E};

        val[0] = (0.006944444444444444*cVal[7])/(dx*dy);
        val[1] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[2] = (0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) + (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[3] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[4] = (-0.006944444444444444*cVal[7])/(dx*dy);
        val[5] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[6] = (0.4444444444444444*cVal[7])/(dx*dy);
        val[7] = (-0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) - (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[8] = (-0.4444444444444444*cVal[7])/(dx*dy);
        val[9] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[10] = (0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) + (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[11] = (-0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) - (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[12] = 1.*cVal[0] - (2.5*cVal[4])/(dx*dx) - (2.5*cVal[5])/(dy*dy) + 1.*cVal[3]*derivativeCoefsZ[1][0] + 1.*cVal[6]*derivativeCoefsZ[2][0];
        val[13] = (0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) + (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[14] = (-0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) - (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[15] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[16] = (-0.4444444444444444*cVal[7])/(dx*dy);
        val[17] = (0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) + (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[18] = (0.4444444444444444*cVal[7])/(dx*dy);
        val[19] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[20] = (-0.006944444444444444*cVal[7])/(dx*dy);
        val[21] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[22] = (-0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) - (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[23] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[24] = (0.006944444444444444*cVal[7])/(dx*dy);
        val[25] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[26] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[27] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[28] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[29] = 1.*cVal[3]*derivativeCoefsZ[1][1] + 1.*cVal[6]*derivativeCoefsZ[2][1];
        val[30] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[31] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[32] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[33] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[34] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[35] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[36] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[37] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[38] = 1.*cVal[3]*derivativeCoefsZ[1][2] + 1.*cVal[6]*derivativeCoefsZ[2][2];
        val[39] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[40] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[41] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[42] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;

        MatSetValuesStencil(jac, 1, &row, 43, col, val, INSERT_VALUES); 

        """.format(eom="{eom}",C="Coefs_2_0_{eom}_", E="_{field}(derivs, mu,Q,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B, c1)")

f = open("HeadersSecondOrder/ExpressionCoefficients_2_0_0.h", "w")

#f.write("""
#        #ifndef __EXP_COEF_2_0_H
#        #define __EXP_COEF_2_0_H
#        """)
f.write(formatExpression_2_0_NoDOF)
f.write("\n")

for n_field in range(num_EOMs):
    f.write(formatExpression_Bdy_DOFOnly.format(n=n_field))

    for n_eom in range(num_EOMs):

        f.write(valExpression_2_0.format(field=n_field, eom=n_eom))
#
#f.write("""
#        #endif
#        """)
f.write("\n")
f.close()

formatExpression_2_1_NoDOF=\
"""
                    col[0].k = k + -2; col[0].j = j + -2; col[0].i = i + 0;
                    col[1].k = k + -2; col[1].j = j + -1; col[1].i = i + 0;
                    col[2].k = k + -2; col[2].j = j + 0; col[2].i = i + -2;
                    col[3].k = k + -2; col[3].j = j + 0; col[3].i = i + -1;
                    col[4].k = k + -2; col[4].j = j + 0; col[4].i = i + 0;
                    col[5].k = k + -2; col[5].j = j + 0; col[5].i = i + 1;
                    col[6].k = k + -2; col[6].j = j + 0; col[6].i = i + 2;
                    col[7].k = k + -2; col[7].j = j + 1; col[7].i = i + 0;
                    col[8].k = k + -2; col[8].j = j + 2; col[8].i = i + 0;
                    col[9].k = k + -1; col[9].j = j + -2; col[9].i = i + 0;
                    col[10].k = k + -1; col[10].j = j + -1; col[10].i = i + 0;
                    col[11].k = k + -1; col[11].j = j + 0; col[11].i = i + -2;
                    col[12].k = k + -1; col[12].j = j + 0; col[12].i = i + -1;
                    col[13].k = k + -1; col[13].j = j + 0; col[13].i = i + 0;
                    col[14].k = k + -1; col[14].j = j + 0; col[14].i = i + 1;
                    col[15].k = k + -1; col[15].j = j + 0; col[15].i = i + 2;
                    col[16].k = k + -1; col[16].j = j + 1; col[16].i = i + 0;
                    col[17].k = k + -1; col[17].j = j + 2; col[17].i = i + 0;
                    col[18].k = k + 0; col[18].j = j + -2; col[18].i = i + -2;
                    col[19].k = k + 0; col[19].j = j + -2; col[19].i = i + -1;
                    col[20].k = k + 0; col[20].j = j + -2; col[20].i = i + 0;
                    col[21].k = k + 0; col[21].j = j + -2; col[21].i = i + 1;
                    col[22].k = k + 0; col[22].j = j + -2; col[22].i = i + 2;
                    col[23].k = k + 0; col[23].j = j + -1; col[23].i = i + -2;
                    col[24].k = k + 0; col[24].j = j + -1; col[24].i = i + -1;
                    col[25].k = k + 0; col[25].j = j + -1; col[25].i = i + 0;
                    col[26].k = k + 0; col[26].j = j + -1; col[26].i = i + 1;
                    col[27].k = k + 0; col[27].j = j + -1; col[27].i = i + 2;
                    col[28].k = k + 0; col[28].j = j + 0; col[28].i = i + -2;
                    col[29].k = k + 0; col[29].j = j + 0; col[29].i = i + -1;
                    col[30].k = k + 0; col[30].j = j + 0; col[30].i = i + 0;
                    col[31].k = k + 0; col[31].j = j + 0; col[31].i = i + 1;
                    col[32].k = k + 0; col[32].j = j + 0; col[32].i = i + 2;
                    col[33].k = k + 0; col[33].j = j + 1; col[33].i = i + -2;
                    col[34].k = k + 0; col[34].j = j + 1; col[34].i = i + -1;
                    col[35].k = k + 0; col[35].j = j + 1; col[35].i = i + 0;
                    col[36].k = k + 0; col[36].j = j + 1; col[36].i = i + 1;
                    col[37].k = k + 0; col[37].j = j + 1; col[37].i = i + 2;
                    col[38].k = k + 0; col[38].j = j + 2; col[38].i = i + -2;
                    col[39].k = k + 0; col[39].j = j + 2; col[39].i = i + -1;
                    col[40].k = k + 0; col[40].j = j + 2; col[40].i = i + 0;
                    col[41].k = k + 0; col[41].j = j + 2; col[41].i = i + 1;
                    col[42].k = k + 0; col[42].j = j + 2; col[42].i = i + 2;
"""



valExpression_2_1=\
        """
        row.c = {eom};

        cVal[0] = {C}0{E};
        cVal[1] = {C}1{E};
        cVal[2] = {C}2{E};
        cVal[3] = {C}3{E};
        cVal[4] = {C}4{E};
        cVal[5] = {C}5{E};
        cVal[6] = {C}6{E};
        cVal[7] = {C}7{E};
        cVal[8] = {C}8{E};
        cVal[9] = {C}9{E};

        val[0] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[1] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[2] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[3] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[4] = 1.*cVal[3]*derivativeCoefsZ[1][0] + 1.*cVal[6]*derivativeCoefsZ[2][0];
        val[5] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[6] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][0])/dx;
        val[7] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[8] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][0])/dy;
        val[9] = (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[10] = (-0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[11] = (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[12] = (-0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[13] = 1.*cVal[3]*derivativeCoefsZ[1][1] + 1.*cVal[6]*derivativeCoefsZ[2][1];
        val[14] = (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[15] = (-0.08333333333333333*cVal[8]*derivativeCoefsZ[1][1])/dx;
        val[16] = (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[17] = (-0.08333333333333333*cVal[9]*derivativeCoefsZ[1][1])/dy;
        val[18] = (0.006944444444444444*cVal[7])/(dx*dy);
        val[19] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[20] = (0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) + (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[21] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[22] = (-0.006944444444444444*cVal[7])/(dx*dy);
        val[23] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[24] = (0.4444444444444444*cVal[7])/(dx*dy);
        val[25] = (-0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) - (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[26] = (-0.4444444444444444*cVal[7])/(dx*dy);
        val[27] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[28] = (0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) + (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[29] = (-0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) - (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[30] = 1.*cVal[0] - (2.5*cVal[4])/(dx*dx) - (2.5*cVal[5])/(dy*dy) + 1.*cVal[3]*derivativeCoefsZ[1][2] + 1.*cVal[6]*derivativeCoefsZ[2][2];
        val[31] = (0.6666666666666666*cVal[1])/dx + (1.3333333333333333*cVal[4])/(dx*dx) + (0.6666666666666666*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[32] = (-0.08333333333333333*cVal[1])/dx - (0.08333333333333333*cVal[4])/(dx*dx) - (0.08333333333333333*cVal[8]*derivativeCoefsZ[1][2])/dx;
        val[33] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[34] = (-0.4444444444444444*cVal[7])/(dx*dy);
        val[35] = (0.6666666666666666*cVal[2])/dy + (1.3333333333333333*cVal[5])/(dy*dy) + (0.6666666666666666*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[36] = (0.4444444444444444*cVal[7])/(dx*dy);
        val[37] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[38] = (-0.006944444444444444*cVal[7])/(dx*dy);
        val[39] = (0.05555555555555555*cVal[7])/(dx*dy);
        val[40] = (-0.08333333333333333*cVal[2])/dy - (0.08333333333333333*cVal[5])/(dy*dy) - (0.08333333333333333*cVal[9]*derivativeCoefsZ[1][2])/dy;
        val[41] = (-0.05555555555555555*cVal[7])/(dx*dy);
        val[42] = (0.006944444444444444*cVal[7])/(dx*dy);
        
        MatSetValuesStencil(jac, 1, &row, 43, col, val, INSERT_VALUES); 

        """.format(eom="{eom}",C="Coefs_2_1_{eom}_", E="_{field}(derivs, mu,Q,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")



f = open("HeadersSecondOrder/ExpressionCoefficients_2_1_0.h", "w")

#f.write("""
#        #ifndef __EXP_COEF_2_1_H
#        #define __EXP_COEF_2_1_H
#        """)
f.write(formatExpression_2_1_NoDOF)
f.write("\n")

for n_field in range(num_EOMs):
    f.write(formatExpression_Bdy_DOFOnly.format(n=n_field))

    for n_eom in range(num_EOMs):

        f.write(valExpression_2_1.format(field=n_field, eom=n_eom))
#f.write("""
#        #endif
#        """)
#
f.write("\n")
f.close()

