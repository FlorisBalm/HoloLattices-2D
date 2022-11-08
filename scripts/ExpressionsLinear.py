num_EOMs=8
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

                    val[0]=cVal[9]/(144.*dy*dz);
                    val[1]=-cVal[9]/(18.*dy*dz);
                    val[2]=cVal[8]/(144.*dx*dz);
                    val[3]=-cVal[8]/(18.*dx*dz);
                    val[4]=cVal[3]/(12.*dz) - cVal[6]/(12.*(dz*dz));
                    val[5]=cVal[8]/(18.*dx*dz);
                    val[6]=-cVal[8]/(144.*dx*dz);
                    val[7]=cVal[9]/(18.*dy*dz);
                    val[8]=-cVal[9]/(144.*dy*dz);
                    val[9]=-cVal[9]/(18.*dy*dz);
                    val[10]=(4*cVal[9])/(9.*dy*dz);
                    val[11]=-cVal[8]/(18.*dx*dz);
                    val[12]=(4*cVal[8])/(9.*dx*dz);
                    val[13]=(-2*cVal[3])/(3.*dz) + (4*cVal[6])/(3.*(dz*dz));
                    val[14]=(-4*cVal[8])/(9.*dx*dz);
                    val[15]=cVal[8]/(18.*dx*dz);
                    val[16]=(-4*cVal[9])/(9.*dy*dz);
                    val[17]=cVal[9]/(18.*dy*dz);
                    val[18]=cVal[7]/(144.*dx*dy);
                    val[19]=-cVal[7]/(18.*dx*dy);
                    val[20]=cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy));
                    val[21]=cVal[7]/(18.*dx*dy);
                    val[22]=-cVal[7]/(144.*dx*dy);
                    val[23]=-cVal[7]/(18.*dx*dy);
                    val[24]=(4*cVal[7])/(9.*dx*dy);
                    val[25]=(-2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy));
                    val[26]=(-4*cVal[7])/(9.*dx*dy);
                    val[27]=cVal[7]/(18.*dx*dy);
                    val[28]=cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx));
                    val[29]=(-2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx));
                    val[30]=cVal[0] - (5*cVal[4])/(2.*(dx*dx)) - (5*cVal[5])/(2.*(dy*dy)) - (5*cVal[6])/(2.*(dz*dz));
                    val[31]=(2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx));
                    val[32]=-cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx));
                    val[33]=cVal[7]/(18.*dx*dy);
                    val[34]=(-4*cVal[7])/(9.*dx*dy);
                    val[35]=(2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy));
                    val[36]=(4*cVal[7])/(9.*dx*dy);
                    val[37]=-cVal[7]/(18.*dx*dy);
                    val[38]=-cVal[7]/(144.*dx*dy);
                    val[39]=cVal[7]/(18.*dx*dy);
                    val[40]=-cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy));
                    val[41]=-cVal[7]/(18.*dx*dy);
                    val[42]=cVal[7]/(144.*dx*dy);
                    val[43]=cVal[9]/(18.*dy*dz);
                    val[44]=(-4*cVal[9])/(9.*dy*dz);
                    val[45]=cVal[8]/(18.*dx*dz);
                    val[46]=(-4*cVal[8])/(9.*dx*dz);
                    val[47]=(2*cVal[3])/(3.*dz) + (4*cVal[6])/(3.*(dz*dz));
                    val[48]=(4*cVal[8])/(9.*dx*dz);
                    val[49]=-cVal[8]/(18.*dx*dz);
                    val[50]=(4*cVal[9])/(9.*dy*dz);
                    val[51]=-cVal[9]/(18.*dy*dz);
                    val[52]=-cVal[9]/(144.*dy*dz);
                    val[53]=cVal[9]/(18.*dy*dz);
                    val[54]=-cVal[8]/(144.*dx*dz);
                    val[55]=cVal[8]/(18.*dx*dz);
                    val[56]=-cVal[3]/(12.*dz) - cVal[6]/(12.*(dz*dz));
                    val[57]=-cVal[8]/(18.*dx*dz);
                    val[58]=cVal[8]/(144.*dx*dz);
                    val[59]=-cVal[9]/(18.*dy*dz);
                    val[60]=cVal[9]/(144.*dy*dz);

                    MatSetValuesStencil(jac, 1, &row, 61, col, val, INSERT_VALUES); 
        """.format(eom="{eom}",C="Coefs_I_{eom}_", E="_{field}(derivs, mu,mu1,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")

f = open("HeadersSecondOrderLinear/ExpressionsCoefficientsInternal.h", "w")
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

                    val[0]=-cVal[9]/(36.*dy*dz);
                    val[1]=(2*cVal[9])/(9.*dy*dz);
                    val[2]=-cVal[8]/(36.*dx*dz);
                    val[3]=(2*cVal[8])/(9.*dx*dz);
                    val[4]=-cVal[3]/(3.*dz) + cVal[6]/(dz*dz);
                    val[5]=(-2*cVal[8])/(9.*dx*dz);
                    val[6]=cVal[8]/(36.*dx*dz);
                    val[7]=(-2*cVal[9])/(9.*dy*dz);
                    val[8]=cVal[9]/(36.*dy*dz);
                    val[9]=cVal[7]/(144.*dx*dy);
                    val[10]=-cVal[7]/(18.*dx*dy);
                    val[11]=cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) - cVal[9]/(24.*dy*dz);
                    val[12]=cVal[7]/(18.*dx*dy);
                    val[13]=-cVal[7]/(144.*dx*dy);
                    val[14]=-cVal[7]/(18.*dx*dy);
                    val[15]=(4*cVal[7])/(9.*dx*dy);
                    val[16]=(-2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) + cVal[9]/(3.*dy*dz);
                    val[17]=(-4*cVal[7])/(9.*dx*dy);
                    val[18]=cVal[7]/(18.*dx*dy);
                    val[19]=cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) - cVal[8]/(24.*dx*dz);
                    val[20]=(-2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) + cVal[8]/(3.*dx*dz);
                    val[21]=cVal[0] - cVal[3]/(2.*dz) - (5*cVal[4])/(2.*(dx*dx)) - (5*cVal[5])/(2.*(dy*dy)) - (2*cVal[6])/(dz*dz);
                    val[22]=(2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) - cVal[8]/(3.*dx*dz);
                    val[23]=-cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) + cVal[8]/(24.*dx*dz);
                    val[24]=cVal[7]/(18.*dx*dy);
                    val[25]=(-4*cVal[7])/(9.*dx*dy);
                    val[26]=(2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) - cVal[9]/(3.*dy*dz);
                    val[27]=(4*cVal[7])/(9.*dx*dy);
                    val[28]=-cVal[7]/(18.*dx*dy);
                    val[29]=-cVal[7]/(144.*dx*dy);
                    val[30]=cVal[7]/(18.*dx*dy);
                    val[31]=-cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) + cVal[9]/(24.*dy*dz);
                    val[32]=-cVal[7]/(18.*dx*dy);
                    val[33]=cVal[7]/(144.*dx*dy);
                    val[34]=cVal[9]/(12.*dy*dz);
                    val[35]=(-2*cVal[9])/(3.*dy*dz);
                    val[36]=cVal[8]/(12.*dx*dz);
                    val[37]=(-2*cVal[8])/(3.*dx*dz);
                    val[38]=cVal[3]/dz + cVal[6]/(dz*dz);
                    val[39]=(2*cVal[8])/(3.*dx*dz);
                    val[40]=-cVal[8]/(12.*dx*dz);
                    val[41]=(2*cVal[9])/(3.*dy*dz);
                    val[42]=-cVal[9]/(12.*dy*dz);
                    val[43]=-cVal[9]/(72.*dy*dz);
                    val[44]=cVal[9]/(9.*dy*dz);
                    val[45]=-cVal[8]/(72.*dx*dz);
                    val[46]=cVal[8]/(9.*dx*dz);
                    val[47]=-cVal[3]/(6.*dz);
                    val[48]=-cVal[8]/(9.*dx*dz);
                    val[49]=cVal[8]/(72.*dx*dz);
                    val[50]=-cVal[9]/(9.*dy*dz);
                    val[51]=cVal[9]/(72.*dy*dz);

                    MatSetValuesStencil(jac, 1, &row, 52, col, val, INSERT_VALUES); 
        """.format(eom="{eom}",C="Coefs_I_{eom}_", E="_{field}(derivs, mu,mu1,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")

f = open("HeadersSecondOrderLinear/ExpressionCoefficients_2_0_1.h", "w")
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

                    val[0]=cVal[9]/(72.*dy*dz);
                    val[1]=-cVal[9]/(9.*dy*dz);
                    val[2]=cVal[8]/(72.*dx*dz);
                    val[3]=-cVal[8]/(9.*dx*dz);
                    val[4]=cVal[3]/(6.*dz);
                    val[5]=cVal[8]/(9.*dx*dz);
                    val[6]=-cVal[8]/(72.*dx*dz);
                    val[7]=cVal[9]/(9.*dy*dz);
                    val[8]=-cVal[9]/(72.*dy*dz);
                    val[9]=-cVal[9]/(12.*dy*dz);
                    val[10]=(2*cVal[9])/(3.*dy*dz);
                    val[11]=-cVal[8]/(12.*dx*dz);
                    val[12]=(2*cVal[8])/(3.*dx*dz);
                    val[13]=-(cVal[3]/dz) + cVal[6]/(dz*dz);
                    val[14]=(-2*cVal[8])/(3.*dx*dz);
                    val[15]=cVal[8]/(12.*dx*dz);
                    val[16]=(-2*cVal[9])/(3.*dy*dz);
                    val[17]=cVal[9]/(12.*dy*dz);
                    val[18]=cVal[7]/(144.*dx*dy);
                    val[19]=-cVal[7]/(18.*dx*dy);
                    val[20]=cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) + cVal[9]/(24.*dy*dz);
                    val[21]=cVal[7]/(18.*dx*dy);
                    val[22]=-cVal[7]/(144.*dx*dy);
                    val[23]=-cVal[7]/(18.*dx*dy);
                    val[24]=(4*cVal[7])/(9.*dx*dy);
                    val[25]=(-2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) - cVal[9]/(3.*dy*dz);
                    val[26]=(-4*cVal[7])/(9.*dx*dy);
                    val[27]=cVal[7]/(18.*dx*dy);
                    val[28]=cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) + cVal[8]/(24.*dx*dz);
                    val[29]=(-2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) - cVal[8]/(3.*dx*dz);
                    val[30]=cVal[0] + cVal[3]/(2.*dz) - (5*cVal[4])/(2.*(dx*dx)) - (5*cVal[5])/(2.*(dy*dy)) - (2*cVal[6])/(dz*dz);
                    val[31]=(2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) + cVal[8]/(3.*dx*dz);
                    val[32]=-cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) - cVal[8]/(24.*dx*dz);
                    val[33]=cVal[7]/(18.*dx*dy);
                    val[34]=(-4*cVal[7])/(9.*dx*dy);
                    val[35]=(2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) + cVal[9]/(3.*dy*dz);
                    val[36]=(4*cVal[7])/(9.*dx*dy);
                    val[37]=-cVal[7]/(18.*dx*dy);
                    val[38]=-cVal[7]/(144.*dx*dy);
                    val[39]=cVal[7]/(18.*dx*dy);
                    val[40]=-cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) - cVal[9]/(24.*dy*dz);
                    val[41]=-cVal[7]/(18.*dx*dy);
                    val[42]=cVal[7]/(144.*dx*dy);
                    val[43]=cVal[9]/(36.*dy*dz);
                    val[44]=(-2*cVal[9])/(9.*dy*dz);
                    val[45]=cVal[8]/(36.*dx*dz);
                    val[46]=(-2*cVal[8])/(9.*dx*dz);
                    val[47]=cVal[3]/(3.*dz) + cVal[6]/(dz*dz);
                    val[48]=(2*cVal[8])/(9.*dx*dz);
                    val[49]=-cVal[8]/(36.*dx*dz);
                    val[50]=(2*cVal[9])/(9.*dy*dz);
                    val[51]=-cVal[9]/(36.*dy*dz);

                    MatSetValuesStencil(jac, 1, &row, 52, col, val, INSERT_VALUES); 
        """.format(eom="{eom}",C="Coefs_I_{eom}_", E="_{field}(derivs, mu,mu1,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")

f = open("HeadersSecondOrderLinear/ExpressionCoefficients_2_1_1.h", "w")
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

        val[0]=cVal[7]/(144.*dx*dy);
        val[1]=-cVal[7]/(18.*dx*dy);
        val[2]=cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) - cVal[9]/(8.*dy*dz);
        val[3]=cVal[7]/(18.*dx*dy);
        val[4]=-cVal[7]/(144.*dx*dy);
        val[5]=-cVal[7]/(18.*dx*dy);
        val[6]=(4*cVal[7])/(9.*dx*dy);
        val[7]=(-2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) + cVal[9]/(dy*dz);
        val[8]=(-4*cVal[7])/(9.*dx*dy);
        val[9]=cVal[7]/(18.*dx*dy);
        val[10]=cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) - cVal[8]/(8.*dx*dz);
        val[11]=(-2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) + cVal[8]/(dx*dz);
        val[12]=cVal[0] - (3*cVal[3])/(2.*dz) - (5*cVal[4])/(2.*(dx*dx)) - (5*cVal[5])/(2.*(dy*dy)) + cVal[6]/(dz*dz);
        val[13]=(2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) - cVal[8]/(dx*dz);
        val[14]=-cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) + cVal[8]/(8.*dx*dz);
        val[15]=cVal[7]/(18.*dx*dy);
        val[16]=(-4*cVal[7])/(9.*dx*dy);
        val[17]=(2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) - cVal[9]/(dy*dz);
        val[18]=(4*cVal[7])/(9.*dx*dy);
        val[19]=-cVal[7]/(18.*dx*dy);
        val[20]=-cVal[7]/(144.*dx*dy);
        val[21]=cVal[7]/(18.*dx*dy);
        val[22]=-cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) + cVal[9]/(8.*dy*dz);
        val[23]=-cVal[7]/(18.*dx*dy);
        val[24]=cVal[7]/(144.*dx*dy);
        val[25]=cVal[9]/(6.*dy*dz);
        val[26]=(-4*cVal[9])/(3.*dy*dz);
        val[27]=cVal[8]/(6.*dx*dz);
        val[28]=(-4*cVal[8])/(3.*dx*dz);
        val[29]=(2*cVal[3])/dz - (2*cVal[6])/(dz*dz);
        val[30]=(4*cVal[8])/(3.*dx*dz);
        val[31]=-cVal[8]/(6.*dx*dz);
        val[32]=(4*cVal[9])/(3.*dy*dz);
        val[33]=-cVal[9]/(6.*dy*dz);
        val[34]=-cVal[9]/(24.*dy*dz);
        val[35]=cVal[9]/(3.*dy*dz);
        val[36]=-cVal[8]/(24.*dx*dz);
        val[37]=cVal[8]/(3.*dx*dz);
        val[38]=-cVal[3]/(2.*dz) + cVal[6]/(dz*dz);
        val[39]=-cVal[8]/(3.*dx*dz);
        val[40]=cVal[8]/(24.*dx*dz);
        val[41]=-cVal[9]/(3.*dy*dz);
        val[42]=cVal[9]/(24.*dy*dz);

        MatSetValuesStencil(jac, 1, &row, 43, col, val, INSERT_VALUES); 

        """.format(eom="{eom}",C="Coefs_2_0_{eom}_", E="_{field}(derivs, mu,mu1,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B, c1)")

f = open("HeadersSecondOrderLinear/ExpressionCoefficients_2_0_0.h", "w")

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

        val[0]=cVal[9]/(24.*dy*dz);
        val[1]=-cVal[9]/(3.*dy*dz);
        val[2]=cVal[8]/(24.*dx*dz);
        val[3]=-cVal[8]/(3.*dx*dz);
        val[4]=cVal[3]/(2.*dz) + cVal[6]/(dz*dz);
        val[5]=cVal[8]/(3.*dx*dz);
        val[6]=-cVal[8]/(24.*dx*dz);
        val[7]=cVal[9]/(3.*dy*dz);
        val[8]=-cVal[9]/(24.*dy*dz);
        val[9]=-cVal[9]/(6.*dy*dz);
        val[10]=(4*cVal[9])/(3.*dy*dz);
        val[11]=-cVal[8]/(6.*dx*dz);
        val[12]=(4*cVal[8])/(3.*dx*dz);
        val[13]=(-2*cVal[3])/dz - (2*cVal[6])/(dz*dz);
        val[14]=(-4*cVal[8])/(3.*dx*dz);
        val[15]=cVal[8]/(6.*dx*dz);
        val[16]=(-4*cVal[9])/(3.*dy*dz);
        val[17]=cVal[9]/(6.*dy*dz);
        val[18]=cVal[7]/(144.*dx*dy);
        val[19]=-cVal[7]/(18.*dx*dy);
        val[20]=cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) + cVal[9]/(8.*dy*dz);
        val[21]=cVal[7]/(18.*dx*dy);
        val[22]=-cVal[7]/(144.*dx*dy);
        val[23]=-cVal[7]/(18.*dx*dy);
        val[24]=(4*cVal[7])/(9.*dx*dy);
        val[25]=(-2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) - cVal[9]/(dy*dz);
        val[26]=(-4*cVal[7])/(9.*dx*dy);
        val[27]=cVal[7]/(18.*dx*dy);
        val[28]=cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) + cVal[8]/(8.*dx*dz);
        val[29]=(-2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) - cVal[8]/(dx*dz);
        val[30]=cVal[0] + (3*cVal[3])/(2.*dz) - (5*cVal[4])/(2.*(dx*dx)) - (5*cVal[5])/(2.*(dy*dy)) + cVal[6]/(dz*dz);
        val[31]=(2*cVal[1])/(3.*dx) + (4*cVal[4])/(3.*(dx*dx)) + cVal[8]/(dx*dz);
        val[32]=-cVal[1]/(12.*dx) - cVal[4]/(12.*(dx*dx)) - cVal[8]/(8.*dx*dz);
        val[33]=cVal[7]/(18.*dx*dy);
        val[34]=(-4*cVal[7])/(9.*dx*dy);
        val[35]=(2*cVal[2])/(3.*dy) + (4*cVal[5])/(3.*(dy*dy)) + cVal[9]/(dy*dz);
        val[36]=(4*cVal[7])/(9.*dx*dy);
        val[37]=-cVal[7]/(18.*dx*dy);
        val[38]=-cVal[7]/(144.*dx*dy);
        val[39]=cVal[7]/(18.*dx*dy);
        val[40]=-cVal[2]/(12.*dy) - cVal[5]/(12.*(dy*dy)) - cVal[9]/(8.*dy*dz);
        val[41]=-cVal[7]/(18.*dx*dy);
        val[42]=cVal[7]/(144.*dx*dy);

        MatSetValuesStencil(jac, 1, &row, 43, col, val, INSERT_VALUES); 

        """.format(eom="{eom}",C="Coefs_2_1_{eom}_", E="_{field}(derivs, mu,mu1,Gx,Gy,ax,ay,x_coord,y_coord,z_coord,nperiodsx,nperiodsy,phasex,phasey, B,c1)")



f = open("HeadersSecondOrderLinear/ExpressionCoefficients_2_1_0.h", "w")

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

f.close()

