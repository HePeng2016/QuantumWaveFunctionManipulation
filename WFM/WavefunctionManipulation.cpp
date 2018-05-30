#include  "WavafunctionManipulation.h"


int  RotateOperator(std::vector<GTFtype>&b2,std::vector< Coordinate >&AtomN,std::vector< double >&MO_Old,std::vector< double >&MO_New,Matrix rotationM)
{
    unsigned int Stack[512];

    for(int i=0;i<AtomN.size();i++)
    {
       double TempX = AtomN[i].X;
       double TempY = AtomN[i].Y;
       double TempZ = AtomN[i].Z;

       AtomN[i].X = Matrix(rotationM,0,0)*TempX+Matrix(rotationM,0,1)*TempY+Matrix(rotationM,0,2)*TempZ;
       AtomN[i].Y = Matrix(rotationM,1,0)*TempX+Matrix(rotationM,1,1)*TempY+Matrix(rotationM,1,2)*TempZ;
       AtomN[i].Z = Matrix(rotationM,2,0)*TempX+Matrix(rotationM,2,1)*TempY+Matrix(rotationM,2,2)*TempZ;
    }


     std::vector<size_t>Index = ordered(b2,true);

      MO_New.resize(b2.size());

      for(int i=0;i<b2.size();i++)
      {
         MO_New[i]=0;
      }

    for(int j=0;j<Index.size();j++)
    {

      int offset=1;
      if(j<(Index.size()-offset))
      {
        while((b2[Index[j+offset]].Center==b2[Index[j]].Center)&&(b2[Index[j]].exps==b2[Index[j+offset]].exps))
        {
            offset++;
        };
      }


        unsigned int OrginalSumSize=-1;
        int TotalSize=0;
        int TypeBase =0;

        for(int i=0;i<offset;i++)
       {
         int index = i+j;


         int ix=typeTx[b2[Index[index]].functype];
         int iy=typeTy[b2[Index[index]].functype];
         int iz=typeTz[b2[Index[index]].functype];

         unsigned int SumSize = ix+iy+iz;

         if(OrginalSumSize!=SumSize)
         {
             TypeBase = TotalSize;
             TotalSize = typeSize[SumSize]+TotalSize;
             OrginalSumSize=SumSize;
         }

          for(int i1=0;i1<=SumSize;i1++)
            for(int i2=0;i2<=SumSize-i1;i2++)
            {
              int i3 = SumSize-i1-i2;
              int E;
              double MOValue=0;
              decipherType(i1,i2,i3,E);
              int Index_New = E-typeSizeBase[SumSize]+j+TypeBase;

             {

                   int I1x,I1y,I1z,I2x,I2y,I2z,I3x,I3y,I3z;

                   for(I1x=0;(I1x<=i1)&&(I1x<=ix);I1x++)
                   {
                       for(I1y=0;(I1y<=iy)&&(I1y<=i1-I1x);I1y++)
                       {
                           I1z = i1-I1x-I1y;
                           if(I1z>iz)
                           {
                              continue;
                           }

                           for(I2x=0;((I2x<=(ix-I1x))&&(I2x<=i2));I2x++)
                           {
                                for(I2y=0;((I2y<=(iy-I1y))&&(I2y<=(i2-I2x)));I2y++)
                                {
                                    I2z = i2-I2x-I2y;
                                    if(I2z>iz-I1z)
                                    {
                                       continue;
                                    }

                                      for(I3x=0;(I3x<=(ix-I1x-I2x))&&(I3x<=i3);I3x++)
                                     {
                                         for(I3y=0;(I3y<=(iy-I1y-I2y))&&(I3y<=i3-I3x);I3y++)
                                         {
                                              I3z = i3-I3x-I3y;

                                              if(I3z>iz-I1z-I2z)
                                              {
                                                 continue;
                                              }


                                               MOValue =  pow(Matrix(rotationM,0,0),I1x);
                                               MOValue = MOValue*pow(Matrix(rotationM,0,1),I2x);
                                               MOValue = MOValue*pow(Matrix(rotationM,0,2),I3x);
                                               MOValue = MOValue*pow(Matrix(rotationM,1,0),I1y);
                                               MOValue = MOValue*pow(Matrix(rotationM,1,1),I2y);
                                               MOValue = MOValue*pow(Matrix(rotationM,1,2),I3y);
                                               MOValue = MOValue*pow(Matrix(rotationM,2,0),I1z);
                                               MOValue = MOValue*pow(Matrix(rotationM,2,1),I2z);
                                               MOValue = MOValue*pow(Matrix(rotationM,2,2),I3z);
                                               CoefficientAttanch(I1x,I2x,I3x,MOValue);
                                               CoefficientAttanch(I1y,I2y,I3y,MOValue);
                                               CoefficientAttanch(I1z,I2z,I3z,MOValue);

                                               MO_New[Index[Index_New]]= MOValue*MO_Old[Index[index]]+MO_New[Index[Index_New]];
                                              printf("%d,%d,%d,%d,%d,%d,%d,%d,%d\n",I1x,I1y,I1z,I2x,I2y,I2z,I3x,I3y,I3z);

                                         }
                                     }
                                }
                           }
                       }
                   }
                   printf("------------\n");
               }


            }

        }
        if(TotalSize!=offset)
        {
          printf(" GFT Type is not compatibility.\n");
          return -1;
        }
        j = j+offset-1;
      }
    Index.resize(0);
    return 1;
}

void  MO_Operator(std::vector<STFtype>&b2,std::vector< std::vector<Operator> >&List,std::vector< Coordinate >&AtomN,unsigned int Operator_Type)
{
     List.resize(b2.size());

    if(Operator_Type== X_P )
    {


        for(int i=0;i<b2.size();i++)
        {
            List[i].resize(2);
            List[i][0].OperatorType=X_phi_;
            List[i][0].OperatorParameter=1.0;
            List[i][1].OperatorType=0;
            List[i][1].OperatorParameter=AtomN[b2[i].Center].X;

        }
    }
    if(Operator_Type== Y_P )
    {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(2);
             List[i][0].OperatorType=Y_phi_;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=0;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Y;

        }
    }

     if(Operator_Type== Z_P)
    {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(2);
             List[i][0].OperatorType=Z_phi_;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=0;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Z;

        }
    }

      if(Operator_Type== XX_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=X_phi_*2;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=X_phi_;
             List[i][1].OperatorParameter=2*AtomN[b2[i].Center].X;
             List[i][2].OperatorType=0;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].X*AtomN[b2[i].Center].X;
        }
     }


      if(Operator_Type== YY_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=Y_phi_*2;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=Y_phi_;
             List[i][1].OperatorParameter=2*AtomN[b2[i].Center].Y;
             List[i][2].OperatorType=0;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].Y*AtomN[b2[i].Center].Y;
        }
     }
     if(Operator_Type== ZZ_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=Z_phi_*2;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=Z_phi_;
             List[i][1].OperatorParameter=2*AtomN[b2[i].Center].Z;
             List[i][2].OperatorType=0;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].Z*AtomN[b2[i].Center].Z;
        }
      }

     if(Operator_Type== XY_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=(X_phi_|Y_phi_);
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=X_phi_;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Y;
             List[i][2].OperatorType=Y_phi_;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].X;
        }
      }

     if(Operator_Type== XZ_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=(X_phi_|Z_phi_);
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=X_phi_;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Z;
             List[i][2].OperatorType=Z_phi_;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].X;
        }
      }

     if(Operator_Type== YZ_P)
     {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=(Y_phi_|Z_phi_);
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=Y_phi_;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Z;
             List[i][2].OperatorType=Z_phi_;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].Y;
        }
      }



}


void  MO_Operator_GTF(std::vector<GTFtype>&b2,std::vector< std::vector<Operator> >&List,std::vector< Coordinate >&AtomN,unsigned int Operator_Type)
{
     List.resize(b2.size());

    if(Operator_Type== X_P )
    {


        for(int i=0;i<b2.size();i++)
        {
            List[i].resize(2);
            List[i][0].OperatorType=X_phi_;
            List[i][0].OperatorParameter=1.0;
            List[i][1].OperatorType=0;
            List[i][1].OperatorParameter=AtomN[b2[i].Center].X;

        }
    }
    if(Operator_Type== Y_P )
    {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(2);
             List[i][0].OperatorType=Y_phi_;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=0;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Y;

        }
    }

     if(Operator_Type== Z_P)
    {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(2);
             List[i][0].OperatorType=Z_phi_;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=0;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Z;

        }
    }

      if(Operator_Type== XX_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=X_phi_*2;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=X_phi_;
             List[i][1].OperatorParameter=2*AtomN[b2[i].Center].X;
             List[i][2].OperatorType=0;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].X*AtomN[b2[i].Center].X;
        }
     }


      if(Operator_Type== YY_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=Y_phi_*2;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=Y_phi_;
             List[i][1].OperatorParameter=2*AtomN[b2[i].Center].Y;
             List[i][2].OperatorType=0;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].Y*AtomN[b2[i].Center].Y;
        }
     }
     if(Operator_Type== ZZ_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=Z_phi_*2;
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=Z_phi_;
             List[i][1].OperatorParameter=2*AtomN[b2[i].Center].Z;
             List[i][2].OperatorType=0;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].Z*AtomN[b2[i].Center].Z;
        }
      }

     if(Operator_Type== XY_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=(X_phi_|Y_phi_);
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=X_phi_;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Y;
             List[i][2].OperatorType=Y_phi_;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].X;
        }
      }

     if(Operator_Type== XZ_P)
      {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=(X_phi_|Z_phi_);
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=X_phi_;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Z;
             List[i][2].OperatorType=Z_phi_;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].X;
        }
      }

     if(Operator_Type== YZ_P)
     {
        for(int i=0;i<b2.size();i++)
        {
             List[i].resize(3);
             List[i][0].OperatorType=(Y_phi_|Z_phi_);
             List[i][0].OperatorParameter=1.0;
             List[i][1].OperatorType=Y_phi_;
             List[i][1].OperatorParameter=AtomN[b2[i].Center].Z;
             List[i][2].OperatorType=Z_phi_;
             List[i][2].OperatorParameter=AtomN[b2[i].Center].Y;
        }
      }



}



double MOsProductOperator_GTF(std::vector< double > MO1,std::vector< double > MO2,std::vector<GTFtype>&b1,std::vector<GTFtype>&b2,std::vector< Coordinate >&AtomN1,std::vector< Coordinate >&AtomN2,std::vector< std::vector<Operator> >&List) //<phi_1|H|phi_2>
{

   double SumProduct=0;



    for(int i=0;i<b1.size();i++)
        for(int j=0;j<b2.size();j++)
    {

       double value=0;


        for(int N3=0;N3<List[i].size();N3++)
        {
           value=value+List[i][N3].OperatorParameter*overlapIntegral((b1[i].functype|List[i][N3].OperatorType<<16),b1[i].exps,b2[j].functype,b2[j].exps,AtomN1[b1[i].Center],AtomN2[b2[j].Center]);
           if(value!=value)
           {
              overlapIntegral((b1[i].functype|List[i][N3].OperatorType<<16),b1[i].exps,b2[j].functype,b2[j].exps,AtomN1[b1[i].Center],AtomN2[b2[j].Center]);
           }

        }
        value = value*MO1[i]*MO2[j];
        SumProduct = SumProduct+value;
    }

     return SumProduct;


}


double normgau(unsigned int itype,double exp)
{
  int ix=typeTx[itype];
  int iy=typeTy[itype];
  int iz=typeTz[itype];

  return pow(((2*exp)/PI),0.75)*sqrt( pow(8*exp,ix+iy+iz)*(ft(ix)*ft(iy)*ft(iz)))/(ft(2*ix)*ft(2*iy)*ft(2*iz));

}

double MOsProductOperator(double *MO1,double * MO2,std::vector<STFtype>&b1,std::vector<STFtype>&b2,std::vector< Coordinate >&AtomN1,std::vector< Coordinate >&AtomN2,std::vector< std::vector<Operator> >&List) //<phi_1|H|phi_2>
{

   double SumProduct=0;



    for(int i=0;i<b1.size();i++)
        for(int j=0;j<b2.size();j++)
    {

       double value=0;
       int N1;
       int N2;

         for(N1=0;N1<b1[i].exps.size();N1++)
           for(N2=0;N2<b2[j].exps.size();N2++)
            {
              for(int N3=0;N3<List[i].size();N3++)
               value=value+List[i][N3].OperatorParameter*b1[i].contracts[N1]*b2[j].contracts[N2]*overlapIntegral((b1[i].functype|List[i][N3].OperatorType<<16),b1[i].exps[N1],b2[j].functype,b2[j].exps[N2],AtomN1[b1[i].Center],AtomN2[b2[j].Center]);
            }
        value = value*MO1[i]*MO2[j];
        SumProduct = SumProduct+value;
    }

     return SumProduct;


}

double MOsProduct(double *MO1,double * MO2,std::vector<STFtype>&b1,std::vector<STFtype>&b2,std::vector< Coordinate >&AtomN1,std::vector< Coordinate >&AtomN2)//<phi_1|phi_2>
{

   double SumProduct=0;



    for(int i=0;i<b1.size();i++)
        for(int j=0;j<b2.size();j++)
    {

       double value=0;
       int N1;
       int N2;
         for(N1=0;N1<b1[i].exps.size();N1++)
           for(N2=0;N2<b2[j].exps.size();N2++)
            {
               value=value+b1[i].contracts[N1]*b2[j].contracts[N2]*overlapIntegral( b1[i].functype,b1[i].exps[N1],b2[j].functype,b2[j].exps[N2],AtomN1[b1[i].Center],AtomN2[b2[j].Center]);


              // double  TestValue=overlapIntegral( b1[i].functype,b1[i].exps[N1],b2[j].functype,b2[j].exps[N2],AtomN1[b1[i].Center],AtomN2[b2[j].Center]);




              /* printf("Integrate[");
               printf("Exp[-");
               printf("%lf",b1[i].exps[N1]);
               printf("*((x1-");
               printf("%lf",AtomN1[b1[i].Center].X);
               printf(")^2+");
               printf("(y1-");
               printf("%lf",AtomN1[b1[i].Center].Y);
               printf(")^2+");
               printf("(z1-");
               printf("%lf",AtomN1[b1[i].Center].Z);
               printf(")^2)-");
               printf("%lf",b2[j].exps[N2]);
               printf("*((x1-");
               printf("%lf",AtomN2[b2[j].Center].X);
               printf(")^2+");
               printf("(y1-");
               printf("%lf",AtomN2[b2[j].Center].Y);
               printf(")^2+");
               printf("(z1-");
               printf("%lf",AtomN2[b2[j].Center].Z);
               printf(")^2)],");
               printf("{x1,-Infinity,Infinity},{y1,-Infinity,Infinity},{z1,-Infinity,Infinity}]\n");

               int  ix1=typeTx[b1[i].functype];
               int  iy1=typeTy[b1[i].functype];
               int  iz1=typeTz[b1[i].functype];
               int  ix2=typeTx[b2[j].functype];
               int  iy2=typeTy[b2[j].functype];
               int  iz2=typeTz[b2[j].functype];

               if((ix1+iy1+iz1+ix2+iy2+iz2)!=0)
               {
                  if(TestValue>0.01)
               {
                   printf("%lf",TestValue);
                   overlapIntegral( b1[i].functype,b1[i].exps[N1],b2[j].functype,b2[j].exps[N2],AtomN1[b1[i].Center],AtomN2[b2[j].Center]);

               }

               }*/


            }

        value = value*MO1[i]*MO2[j];
        value = value;
        SumProduct = SumProduct+value;
    }

     return SumProduct;

}


double MOsProductGTF( vector < double > MO1,vector < double > MO2,std::vector<GTFtype>&b1,std::vector<GTFtype>&b2,std::vector< Coordinate >&AtomN1,std::vector< Coordinate >&AtomN2)//<phi_1|phi_2>
{

   double SumProduct=0;



    for(int i=0;i<b1.size();i++)
        for(int j=0;j<b2.size();j++)
    {

       double value=0;
       int N1;
       int N2;

        value=value+overlapIntegral( b1[i].functype,b1[i].exps,b2[j].functype,b2[j].exps,AtomN1[b1[i].Center],AtomN2[b2[j].Center]);

        value = value*MO1[i]*MO2[j];
        value = value;
        SumProduct = SumProduct+value;
    }

     return SumProduct;

}






double overlapIntegral(unsigned int itype1,double exp1,unsigned int itype2,double exp2,Coordinate A ,Coordinate B)
{

   double rp = exp1+exp2;
   double np = (exp1*exp2)/rp;

  int OperatorType = itype1;
  OperatorType=(OperatorType>>16);
  itype1 = itype1&0xFFFF;

int  ix1=typeTx[itype1];
int  iy1=typeTy[itype1];
int  iz1=typeTz[itype1];
int  ix2=typeTx[itype2];
int  iy2=typeTy[itype2];
int  iz2=typeTz[itype2];

   ix1=ix1+OperatorType&0x1F;
   iy1=iy1+(OperatorType>>6)&0x1F;
   iz1=iz1+(OperatorType>>11)&0x1F;






   double Sx=0;
   double Sy=0;
   double Sz=0;

   for(int i1=0;i1<=(ix1/2);i1++)
   {
       for(int i2=0;i2<=(ix2/2);i2++)
       {

           for(int o=0;o<=((ix1+ix2-2*(i1+i2))/2);o++)
           {
               double value;

                value=(ft(ix1+ix2-2*(i1+i2))*pow(exp1,ix2-i1-2*i2-o)*pow(exp2,ix1-2*i1-i2-o))/(pow(4,i1+i2+o)*ft(i1)*ft(i2)*ft(o));
                if(o%2==1)
                {
                  value=-value;
                }
                value=(pow(rp,2*(i1+i2)+o)*pow((A.X-B.X),ix1+ix2-2*(i1+i2)-2*o)/(ft(ix1-2*i1)*ft(ix2-2*i2)*ft(ix1+ix2-2*(i1+i2)-2*o)))*value;
                Sx=Sx+value;
           }
       }
   }


   if(ix1%2==1)
   {
       Sx=-Sx;
   }

   Sx=ft(ix1)*ft(ix2)*Sx;
   Sx=Sx/pow(rp,ix1+ix2);

  for(int i1=0;i1<=(iy1/2);i1++)
   {
       for(int i2=0;i2<=(iy2/2);i2++)
       {
           for(int o=0;o<=(iy1+iy2-2*(i1+i2))/2;o++)
           {
               double value;

                value=(ft(iy1+iy2-2*(i1+i2))*pow(exp1,iy2-i1-2*i2-o)*pow(exp2,iy1-2*i1-i2-o))/(pow(4,i1+i2+o)*ft(i1)*ft(i2)*ft(o));
                if(o%2==1)
                {
                  value=-value;
                }
                value=pow(rp,2*(i1+i2)+o)*pow((A.Y-B.Y),iy1+iy2-2*(i1+i2)-2*o)/(ft(iy1-2*i1)*ft(iy2-2*i2)*ft(iy1+iy2-2*(i1+i2)-2*o))*value;
                Sy=Sy+value;
           }
       }
   }

   if(iy1%2==1)
   {
      Sy=-Sy;
   }

   Sy=ft(iy1)*ft(iy2)*Sy;
   Sy=Sy/pow(rp,iy1+iy2);

   for(int i1=0;i1<=(iz1/2);i1++)
   {
       for(int i2=0;i2<=(iz2/2);i2++)
       {
           for(int o=0;o<=(iz1+iz2-2*(i1+i2))/2;o++)
           {
               double value;

                value=(ft(iz1+iz2-2*(i1+i2))*pow(exp1,iz2-i1-2*i2-o)*pow(exp2,iz1-2*i1-i2-o))/(pow(4,i1+i2+o)*ft(i1)*ft(i2)*ft(o));
                if(o%2==1)
                {
                  value=-value;
                }
                value=pow(rp,2*(i1+i2)+o)*pow((A.Z-B.Z),iz1+iz2-2*(i1+i2)-2*o)/(ft(iz1-2*i1)*ft(iz2-2*i2)*ft(iz1+iz2-2*(i1+i2)-2*o))*value;
                Sz=Sz+value;
           }
       }
   }

    if(iz1%2==1)
   {
      Sz=-Sz;
   }


   Sz=ft(iz1)*ft(iz2)*Sz;
   Sz=Sz/pow(rp,iz1+iz2);



return  pow(PI/rp,1.5)*exp(-np*(pow((A.X-B.X),2)+ pow((A.Y-B.Y),2)+pow((A.Z-B.Z),2)))*Sx*Sy*Sz;

}
void read31file(char * fileName,std::vector< Coordinate > &AtomN,std::vector<STFtype>&b)
{

    FILE * file31=fopen(fileName,"r");
    char Buffer[1024];
    int  ncenter;
    int  nshell;
    int nprimshell;
    int nprims;
    int nmo;
    int iGTF=0;
    int ibasis=0;


    std::vector<int> shellcon;
    std::vector<int> shell2atom;
    std::vector<int> shellnumbas;
    std::vector<int> shell2prmshl;
    std::vector<double> prmshlexp;
    std::vector<double> cs;
    std::vector<double> cp;
    std::vector<double> cd;
    std::vector<double> cf;
    std::vector<double> cg;
    std::vector<int> bastype;


    fscanf(file31,"%[^\n]",Buffer);
    fscanf(file31,"\n",Buffer);
    fscanf(file31,"%[^\n]",Buffer);
    fscanf(file31,"\n",Buffer);
    fscanf(file31,"%[^\n]",Buffer);
    fscanf(file31,"\n",Buffer);


    fscanf(file31,"%d",&ncenter);
    fscanf(file31,"%d",&nshell);
    fscanf(file31,"%d",&nprimshell);
    fscanf(file31,"\n",Buffer);
    fscanf(file31,"%[^\n]",Buffer);
    fscanf(file31,"\n",Buffer);



    AtomN.resize(ncenter);
    shellcon.resize(nshell);
    shell2atom.resize(nshell);
    shellnumbas.resize(nshell);
    shell2prmshl.resize(nshell);
    prmshlexp.resize(nprimshell);
    cs.resize(nprimshell);
    cp.resize(nprimshell);
    cd.resize(nprimshell);
    cf.resize(nprimshell);
    cg.resize(nprimshell);
    bastype.resize(nshell*15);


    for(int i=0;i<ncenter;i++)
    {

        fscanf(file31,"%d",&AtomN[i].Index);
        fscanf(file31,"%lf",&AtomN[i].X);
        fscanf(file31,"%lf",&AtomN[i].Y);
        fscanf(file31,"%lf",&AtomN[i].Z);
        fscanf(file31,"\n",Buffer);

    }

    fscanf(file31,"%[^\n]",Buffer);
    fscanf(file31,"\n",Buffer);

    int bastypeI=0;

    for(int i=0;i<nshell;i++)
    {
       fscanf(file31,"%d",&shell2atom[i]);
       fscanf(file31,"%d",&shellnumbas[i]);
       fscanf(file31,"%d",&shell2prmshl[i]);
       fscanf(file31,"%d",&shellcon[i]);

       fscanf(file31,"\n",Buffer);

       for(int j=0;j<shellnumbas[i];j++)
       {
          fscanf(file31,"%d",&bastype[bastypeI+j]);
       }
       bastypeI=bastypeI+shellnumbas[i];

       fscanf(file31,"\n",Buffer);
    }
      bastype.resize(bastypeI);

       fscanf(file31,"%[^\n]",Buffer);
       printf("%s",Buffer);
       fscanf(file31,"\n",Buffer);

       for(int i=0;i<nprimshell;i++)
       {
         fscanf(file31,"%lf",&prmshlexp[i]);
         fscanf(file31,"\n",Buffer);
       }

        for(int i=0;i<nprimshell;i++)
       {
         fscanf(file31,"%lf",&cs[i]);
         fscanf(file31,"\n",Buffer);
       }

        for(int i=0;i<nprimshell;i++)
       {
         fscanf(file31,"%lf",&cp[i]);
         fscanf(file31,"\n",Buffer);
       }

        for(int i=0;i<nprimshell;i++)
       {
         fscanf(file31,"%lf",&cf[i]);
         fscanf(file31,"\n",Buffer);
       }


       fscanf(file31,"\n",Buffer);
       if(!feof(file31))
       {
         for(int i=0;i<nprimshell;i++)
        {
          fscanf(file31,"%lf",&cg[i]);
          fscanf(file31,"\n",Buffer);
        }
       }

       nprims=0;
       for(int i=0;i<nshell;i++)
       {
           nprims= nprims+shellcon[i]*shellnumbas[i];
       }

       int b_size =0;

       for(int i=0;i<nshell;i++)
       {
           b_size= b_size+shellnumbas[i];
       }



       b.resize(b_size);

       double valnorm31;
       double valnormnew;

       for(int i=0;i<nshell;i++)
       {
          for(int j=0;j<shellnumbas[i];j++)
          {
              b[i+j].Center = shell2atom[i]-1;
              b[i+j].exps.resize(shellcon[i]);
              b[i+j].contracts.resize(shellcon[i]);


              switch(bastype[i+j])
              {
                  case(1):
                   b[i+j].functype=S_T;
                   break;
                  case(101):
                   b[i+j].functype=X_T;
                   break;
                  case(102):
                   b[i+j].functype=Y_T;
                   break;
                  case(103):
                   b[i+j].functype=Z_T;
                   break;
                  case(201):
                   b[i+j].functype=XX_T;
                   break;
                  case(202):
                   b[i+j].functype=XY_T;
                   break;
                  case(203):
                   b[i+j].functype=XZ_T;
                   break;
                  case(204):
                   b[i+j].functype=YY_T;
                   break;
                  case(205):
                   b[i+j].functype=YZ_T;
                   break;
                  case(206):
                   b[i+j].functype=ZZ_T;
                   break;
                  case(301):
                   b[i+j].functype=XXX_T;
                   break;
                  case(302):
                   b[i+j].functype=XXY_T;
                   break;
                  case(303):
                   b[i+j].functype=XXZ_T;
                   break;
                  case(304):
                   b[i+j].functype=XYY_T;
                   break;
                  case(305):
                   b[i+j].functype=XYZ_T;
                   break;
                  case(306):
                   b[i+j].functype=XZZ_T;
                   break;
                  case(307):
                    b[i+j].functype=YYY_T;
                    break;
                  case(308):
                    b[i+j].functype=YYZ_T;
                    break;
                  case(309):
                    b[i+j].functype=YZZ_T;
                    break;
                  case(310):
                    b[i+j].functype=YZZ_T;
                    break;

              }
              for(int n=0;n<shellcon[i];n++)
              {
                 b[i+j].exps[n]=prmshlexp[shell2prmshl[i]+n-1];


               if(bastype[i]/100==STYPE)
               {
                  b[i+j].contracts[n]=cs[shell2prmshl[i]+n-1];
               }

               if(bastype[i]/100==PTYPE)
               {
                  b[i+j].contracts[n]=cp[shell2prmshl[i]+n-1];
               }
              if(bastype[i]/100==DTYPE)
               {
                   b[i+j].contracts[n]=cd[shell2prmshl[i]+n-1];
                  // b[ibasis].contracts[n]=b[ibasis].contracts[n]*normgau(b[ibasis].functype,b[ibasis].exps[n]);

                   valnorm31=normgau(4,b[i+j].exps[n]);
                   valnormnew=normgau(7,b[i+j].exps[n]);
                   b[i+j].contracts[n]=(b[i+j].contracts[n]/valnorm31)*valnormnew;

               }
                if(bastype[i]/100==FTYPE)
                {
                    b[i+j].contracts[n]=cf[shell2prmshl[i]+n-1];
                  //  b[ibasis].contracts[n]=b[ibasis].contracts[n]*normgau(b[ibasis].functype,b[ibasis].exps[n]);
                   	if (bastype[i]!=301&&bastype[i]!=307&&bastype[i]!=310)
                    {
                      valnorm31=normgau(10,b[i+j].exps[n]);
                      if (bastype[i]==302||bastype[i]==303||bastype[i]==304||bastype[i]==306||bastype[i]==308||bastype[i]==309)
					  {
                         valnormnew=normgau(13,b[i+j].exps[n]);
					  }else if(bastype[i]==305)
					  {
                         valnormnew=normgau(19,b[i+j].exps[n]);
					  }
                    }
                    b[i+j].contracts[n]=(b[i+j].contracts[n]/valnorm31)*valnormnew;
                }
              }
          }
       }
   fclose(file31);
}

void readfWFN(char * fileName,std::vector< Coordinate > &AtomN,std::vector<GTFtype>&b,std::vector< std::vector<double> >&MOs)
{
      char Buffer[1024];
      int MolOrbitalN;
      int PrimitiveN;
      int nCenter;


      FILE * fileWFN=fopen(fileName,"r");

      fscanf(fileWFN,"%[^\n]",Buffer);
      fscanf(fileWFN,"\n",Buffer);
      fscanf(fileWFN,"%s",Buffer);
      fscanf(fileWFN,"%d",&MolOrbitalN);
      fscanf(fileWFN," ",Buffer);
      fscanf(fileWFN,"MOL ORBITALS",Buffer);
      fscanf(fileWFN,"%d",&PrimitiveN);
      fscanf(fileWFN," ",Buffer);
      fscanf(fileWFN,"PRIMITIVES",Buffer);
      fscanf(fileWFN,"%d",&nCenter);
      fscanf(fileWFN," ",Buffer);
      fscanf(fileWFN,"NUCLEI",Buffer);
      fscanf(fileWFN,"\n",Buffer);

      AtomN.resize(nCenter);

      for(int i=0;i<nCenter;i++)
      {
          double charge;

          fscanf(fileWFN,"\n",Buffer);
          fscanf(fileWFN,"%s",Buffer);
          fscanf(fileWFN,"%s",Buffer);
          fscanf(fileWFN,"%s",Buffer);
          fscanf(fileWFN,"%s",Buffer);
          fscanf(fileWFN,"%lf",&AtomN[i].X);
          fscanf(fileWFN,"%lf",&AtomN[i].Y);
          fscanf(fileWFN,"%lf",&AtomN[i].Z);

          fscanf(fileWFN,"%s",Buffer);
          fscanf(fileWFN," ",Buffer);
          fscanf(fileWFN,"=",Buffer);
          fscanf(fileWFN,"  ",Buffer);
          fscanf(fileWFN,"%lf",&charge);
          AtomN[i].Index=(int)charge;
          fscanf(fileWFN,"\n",Buffer);


      }
      b.resize(PrimitiveN);
      fscanf(fileWFN,"\n",Buffer);
      int oldvalue;
      double floatvalue;

      for(int i=0;i<PrimitiveN;i++)
      {
          fscanf(fileWFN,"\n",Buffer);
          fscanf(fileWFN,"CENTRE",Buffer);
          fscanf(fileWFN," ",Buffer);
          fscanf(fileWFN,"ASSIGNMENTS",Buffer);
          fscanf(fileWFN," ",Buffer);
          fscanf(fileWFN,"%d",&oldvalue);
          b[i].Center=oldvalue-1;
          fscanf(fileWFN,"\n",Buffer);
      }
      for(int i=0;i<PrimitiveN;i++)
      {
          fscanf(fileWFN," ",Buffer);
          fscanf(fileWFN,"TYPE ASSIGNMENTS",Buffer);
          fscanf(fileWFN,"%d",&oldvalue);
          b[i].functype=oldvalue-1;
          fscanf(fileWFN,"\n",Buffer);
      }
      for(int i=0;i<PrimitiveN;i++)
      {   fscanf(fileWFN,"\n",Buffer);
          fscanf(fileWFN," ",Buffer);
          fscanf(fileWFN,"EXPONENTS",Buffer);

          fscanf(fileWFN,"  ",Buffer);
          fscanf(fileWFN,"%lf",&floatvalue);
          b[i].exps=floatvalue;
          fscanf(fileWFN,"\n",Buffer);
      }

      MOs.resize(MolOrbitalN);

      for(int i=0;i<MolOrbitalN;i++)
      {
          fscanf(fileWFN,"%[^\n]",Buffer);
          fscanf(fileWFN,"\n",Buffer);

          MOs[i].resize(PrimitiveN);

          for(int j=0;j<PrimitiveN;j++)
          {
              fscanf(fileWFN,"%lf",&MOs[i][j]);
              fscanf(fileWFN,"\n",Buffer);
          }
      }


}


void InitInvariantContext(gsl_Type * Context)
{

    Context->m;
    Context->evec;
    Context->w;
    Context->eval;
    Context->perm;
      Context->m = gsl_matrix_alloc (3, 3);
     Context->evec = gsl_matrix_alloc(3,3);
   Context->w = gsl_eigen_symmv_alloc(3*4);
  Context->eval = gsl_vector_alloc(3);
  Context->perm = gsl_permutation_alloc(3);
  Context->rotationM.col_size=3;
  Context->rotationM.row_size=3;
  Matrix_Init(Context->rotationM);

}
void FreeInvariantContext(gsl_Type * Context)
{
  gsl_eigen_symmv_free(Context->w);
  gsl_vector_free(Context->eval);
  gsl_matrix_free(Context->evec);
  gsl_matrix_free(Context->m);
  gsl_permutation_free(Context->perm);
  Matrix_Free(Context->rotationM);

}


 void  InvariantDirection(std::vector<GTFtype>bs,std::vector<double> MO_i,std::vector< Coordinate > AtomN_i,std::vector<double> MO_o,std::vector< Coordinate > AtomN_o,gsl_Type * Context)
 {


      std::vector< std::vector<Operator> >Operator_List;
      MO_Operator_GTF(bs,Operator_List,AtomN_i,X_P);
      double X_bar = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_i,AtomN_i,Operator_List);
      Operator_List.resize(0);
      MO_Operator_GTF(bs,Operator_List,AtomN_i,Y_P);
      double Y_bar = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_i,AtomN_i,Operator_List);
      Operator_List.resize(0);
      MO_Operator_GTF(bs,Operator_List,AtomN_i,Z_P);
      double Z_bar = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_i,AtomN_i,Operator_List);
      Operator_List.resize(0);

      AtomN_o.resize(AtomN_i.size());

      for(int i=0;i<AtomN_i.size();i++)
     {
       AtomN_o[i].X=AtomN_i[i].X-X_bar;
       AtomN_o[i].Y=AtomN_i[i].Y-Y_bar;
       AtomN_o[i].Z=AtomN_i[i].Z-Z_bar;
     }

        MO_Operator_GTF(bs,Operator_List,AtomN_o,XX_P);
        double XX = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_o,AtomN_o,Operator_List);
        Operator_List.resize(0);
        MO_Operator_GTF(bs,Operator_List,AtomN_o,XY_P);
        double XY = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_o,AtomN_o,Operator_List);
        Operator_List.resize(0);
        MO_Operator_GTF(bs,Operator_List,AtomN_o,YY_P);
        double YY = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_o,AtomN_o,Operator_List);
        Operator_List.resize(0);
        MO_Operator_GTF(bs,Operator_List,AtomN_o,YZ_P);
        double YZ = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_o,AtomN_o,Operator_List);
        Operator_List.resize(0);
        MO_Operator_GTF(bs,Operator_List,AtomN_o,ZZ_P);
        double ZZ = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_o,AtomN_o,Operator_List);
        Operator_List.resize(0);
        MO_Operator_GTF(bs,Operator_List,AtomN_o,XZ_P);
        double XZ = MOsProductOperator_GTF(MO_i,MO_i,bs,bs,AtomN_o,AtomN_o,Operator_List);

         gsl_matrix_set(Context->m,0,0,XX);
         gsl_matrix_set(Context->m,0,1,XY);
         gsl_matrix_set(Context->m,0,2,XZ);
         gsl_matrix_set(Context->m,1,0,XY);
         gsl_matrix_set(Context->m,1,1,YY);
         gsl_matrix_set(Context->m,1,2,YZ);
         gsl_matrix_set(Context->m,2,0,XZ);
         gsl_matrix_set(Context->m,2,1,YZ);
         gsl_matrix_set(Context->m,2,2,ZZ);

         gsl_eigen_symmv(Context->m,Context->eval,Context->evec,Context->w);
         gsl_sort_vector_index(Context->perm,Context->eval);
         gsl_permutation_reverse(Context->perm);
         RotateOperator(bs,AtomN_o,MO_i,MO_o,Context->rotationM);

 };






