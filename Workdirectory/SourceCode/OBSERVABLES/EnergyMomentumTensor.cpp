#ifndef __ENERGYDENSITYMAP__CPP__
#define __ENERGYDENSITYMAP__CPP__

namespace Observables {
        
    namespace EnergyMomentumTensor{
		
		//COMPUTES COMPONENT OF ENERGY MOMENTUM TENSOR BY AVERAGING OVER LONGITUDINAL PLANE FOR EVENT MONITOR //

        void CreateMapBfCollision(std::string fname,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
            
            //CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
            DOUBLE cB[Lattice::Dimension];  DOUBLE cE[Lattice::Dimension];
            
            //SET CONSTANTS
            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                cB[mu]=SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                cE[mu]=SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
            }
			
            //BUFFERS FOR AVERAGE FIELD STRENGTH
            SET_AVG_FIELD_STRENGTH_BUFFERS();
            
            //CREATE OUTPUT STREAM
            std::ofstream OutStream;
            
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            OutStream.open(OutputFile.c_str());
			
			OutStream.precision(OUTPUT_PRECISION);
			
			for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
					
					DOUBLE T00Left=0.0;		DOUBLE T00Right=0.0;

					for(INT z=zLow; z<=zHigh; z++){
						
						DOUBLE EE0Sqr=0.0;   DOUBLE EE1Sqr=0.0;    DOUBLE EE2Sqr=0.0;
						DOUBLE BB0Sqr=0.0;   DOUBLE BB1Sqr=0.0;    DOUBLE BB2Sqr=0.0;
						
						//CALCULATES E'S FOR SLICE t+dt/2 WHILE B'S FOR t+dt
						COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z,U);
						
						//ADD SQUARED FIELD STRENGTH TO LOCAL OBSERVABLES
						BB0Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uyz);
						BB1Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uzx);
						BB2Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uxy);

						for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
		
							EE0Sqr+=SQR(E->Get(x,y,z,0,a)[0]);
							EE1Sqr+=SQR(E->Get(x,y,z,1,a)[0]);
							EE2Sqr+=SQR(E->Get(x,y,z,2,a)[0]);
							
						}
					
						if(z<=(zHigh-1)/2){
							T00Left+=DOUBLE(0.5)*(cE[0]*EE0Sqr+cE[1]*EE1Sqr+cE[2]*EE2Sqr+cB[0]*BB0Sqr+cB[1]*BB1Sqr+cB[2]*BB2Sqr)/Lattice::N[2];

						}
						else{
							T00Right+=DOUBLE(0.5)*(cE[0]*EE0Sqr+cE[1]*EE1Sqr+cE[2]*EE2Sqr+cB[0]*BB0Sqr+cB[1]*BB1Sqr+cB[2]*BB2Sqr)/Lattice::N[2];
						}
					}
					
					OutStream << x << "  " << y << "  " << T00Left << "  "<< T00Right << std::endl;

				}
			}
			
		 //CLOSE OUTPUT STREAM
		 OutStream.close();
		}
		
        //CREATE ENERGY DENSITY MAP
        void CreateMapBfCollision(std::string fname,GaugeLinks *U,ElectricFields *E){
            
            CreateMapBfCollision(fname,0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
            
        }
        
        
        void CreateMapBfCollision(std::string fname){
            
            CreateMapBfCollision(fname,SUNcGaugeLinks::U,SUNcElectricFields::E);
            
        }
		
		void CreateMap(std::string fname, INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,GaugeLinks *UOld,ElectricFields *E){
			
			//CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
			DOUBLE qB[Lattice::Dimension];  DOUBLE qE[Lattice::Dimension];	DOUBLE OJcoef[Lattice::Dimension];
			DOUBLE NormFactor=(1.0/(Lattice::N[0]*Lattice::N[1]));
			//SET CONSTANTS
			for(int mu=0;mu<Lattice::Dimension;mu++){
				OJcoef[mu]=(1.0/(2.0*U->a[mu]*U->aCube));
				qB[mu]=(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
				qE[mu]=(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
			}
			
			
			//CREATE OUTPUT STREAM
			std::ofstream OutStream;

			OutStream.precision(OUTPUT_PRECISION);

			std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
			
			OutStream.open(OutputFile.c_str());

			// Refresh the boundary planes once before the parallel profile reads.
			UOld->SynchronizeGhostCells();
			U->SynchronizeGhostCells();
						
			//AVERAGE THE EMT OVER EACH TRANSVERSE PLANE
			for(INT z=zLow;z<=zHigh;z++){
				//DIAGONAL COMPONENTS OF THE ENERGY MOMENTUM TENSOR
				DOUBLE T00=0.0,TXX=0.0,TYY=0.0,TZZ=0.0;
				DOUBLE T0X=0.0,T0Y=0.0,T0Z=0.0;//,TXY=0.0,TYZ=0.0,TZX=0.0;
				DOUBLE T0Xmod=0.0,T0Ymod=0.0,T0Zmod=0.0;
				DOUBLE TXXmod=0.0, TYYmod=0.0, TZZmod=0.0; 
				 
				

			#pragma omp parallel
				{

				//BUFFERS FOR AVERAGE FIELD STRENGTH
				SET_AVG_FIELD_STRENGTH_BUFFERS(); 
					
			#pragma omp for reduction( + : T00,TXX,TYY,TZZ,T0X,T0Y,T0Z,T0Xmod,T0Ymod,T0Zmod,TXXmod,TYYmod,TZZmod) //TXY,TYZ, TZX removed

				for(INT y=yLow;y<=yHigh;y++){
					for(INT x=xLow;x<=xHigh;x++){
						
						//CALCULATES E'S AND B'S FOR SLICE t+dt/2
						COMPUTE_AVG_FIELD_STRENGTH(x,y,z,UOld);

						//COMPUTE DIAGONAL COMPONENTS OF ENERGY MOMENTUM TENSOR
						T00+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc+SQR(qE[2])*E2SqrLoc)+(SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc+SQR(qB[2])*B2SqrLoc)))/(Lattice::N[0]*Lattice::N[1]);
						
						TXX+=(DOUBLE(0.5)*((SQR(qE[1])*E1SqrLoc+SQR(qE[2])*E2SqrLoc)-(SQR(qE[0])*E0SqrLoc)))*(NormFactor);
						TXX+=(DOUBLE(0.5)*((SQR(qB[1])*B1SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qB[0])*B0SqrLoc)))*(NormFactor);

						//TYY+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qE[1])*E1SqrLoc+SQR(qB[1])*B1SqrLoc)))/(Lattice::N[0]*Lattice::N[1]);
						TYY+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[2])*E2SqrLoc)-(SQR(qE[1])*E1SqrLoc)))*(NormFactor);
						TYY+=(DOUBLE(0.5)*((SQR(qB[0])*B0SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qB[1])*B1SqrLoc)))*(NormFactor);

						//TZZ+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc)-(SQR(qE[2])*E2SqrLoc+SQR(qB[2])*B2SqrLoc)))/(Lattice::N[0]*Lattice::N[1]);
						TZZ+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc)-(SQR(qE[2])*E2SqrLoc)))*(NormFactor);
						TZZ+=(DOUBLE(0.5)*((SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc)-(SQR(qB[2])*B2SqrLoc)))*(NormFactor);
						

						//MAGNETIC PART OF TXXmod
						TXXmod+=(SQR(qB[2]) * (ReTrIDMinusUxympx+ReTrIDMinusUxym+ReTrIDMinusUxympxmy+ReTrIDMinusUxymmy) * NormFactor);
						TXXmod+=(SQR(qB[1]) * (ReTrIDMinusUzxmpx+ReTrIDMinusUzxm+ReTrIDMinusUzxmpxmz+ReTrIDMinusUzxmmz) * NormFactor);
						//T00B part of TXXmod
						TXXmod+=-(2.0 * (SQR(qB[2])) * B2SqrLocmpx * NormFactor);
						TXXmod+=-(2.0 * (SQR(qB[1])) * B1SqrLocmpx * NormFactor);
						TXXmod+=-(2.0 * (SQR(qB[0])) * B0SqrLocmpx * NormFactor);


						//MAGNETIC PART OF TYYmod
						TYYmod+=(SQR(qB[2]) * (ReTrIDMinusUxympy+ReTrIDMinusUxym+ReTrIDMinusUxympymx+ReTrIDMinusUxymmx) * NormFactor);
						TYYmod+=(SQR(qB[0]) * (ReTrIDMinusUyzmpy+ReTrIDMinusUyzm+ReTrIDMinusUyzmpymz+ReTrIDMinusUyzmmz) * NormFactor);
						//T00 part of TYYmod
						TYYmod+=-(2.0 * (SQR(qB[2])) * B2SqrLocmpy * NormFactor);
						TYYmod+=-(2.0 * (SQR(qB[1])) * B1SqrLocmpy * NormFactor);
						TYYmod+=-(2.0 * (SQR(qB[0])) * B0SqrLocmpy * NormFactor);


						//MAGNETIC part of TZZmod
						TZZmod+=((SQR(qB[1])) * (ReTrIDMinusUzxmpz+ReTrIDMinusUzxm+ReTrIDMinusUzxmpzmx+ReTrIDMinusUzxmmx) * NormFactor);
						TZZmod+=((SQR(qB[0])) * (ReTrIDMinusUyzmpz+ReTrIDMinusUyzm+ReTrIDMinusUyzmpzmy+ReTrIDMinusUyzmmy) * NormFactor);
						//T00 part of TZZmod
						TZZmod+=-(2.0 * (SQR(qB[2])) * B2SqrLocmpz * NormFactor);
						TZZmod+=-(2.0 * (SQR(qB[1])) * B1SqrLocmpz * NormFactor);
						TZZmod+=-(2.0 * (SQR(qB[0])) * B0SqrLocmpz * NormFactor);


						//COMPUTE POYNTING VECTOR
						T0X+=(qE[1]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B2Loc)-qE[2]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B1Loc))*(NormFactor);
						T0Y+=(qE[2]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B0Loc)-qE[0]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B2Loc))*(NormFactor);
						T0Z+=(qE[0]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B1Loc)-qE[1]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B0Loc))*(NormFactor);

						//COMPUTE SPATIAL STRESS
						// TXY+=(qE[0]*qE[1]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,E1Loc)+qB[0]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(B0Loc,B1Loc))*(NormFactor);
						// TYZ+=(qE[1]*qE[2]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,E2Loc)+qB[1]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(B1Loc,B2Loc))*(NormFactor);
						// TZX+=(qE[2]*qE[0]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,E0Loc)+qB[2]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(B2Loc,B0Loc))*(NormFactor);

						//for-loop for colours
						for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

							//TXXmod Electric part//
							TXXmod+=((1.0 / (4.0 * SQR(U->a[0]) * SQR(U->a[2]))) * (SQR(E->Get(x+1,y,z,1,a)[0]) + SQR(E->Get(x+1,y-1,z,1,a)[0])) * NormFactor);
							TXXmod+=((1.0 / (4.0 * SQR(U->a[0]) * SQR(U->a[1]))) * (SQR(E->Get(x+1,y,z,2,a)[0]) + SQR(E->Get(x+1,y,z-1,2,a)[0])) * NormFactor);
							TXXmod+=-((1.0 / (SQR(U->a[1]) * SQR(U->a[2]))) * 0.5 * (E0Up[a] * E->Get(x+1,y,z,0,a)[0]) * NormFactor);

							//TYYmod Electric part//
							TYYmod+=((1.0 / (4.0 * SQR(U->a[1]) * SQR(U->a[2]))) * (SQR(E->Get(x,y+1,z,0,a)[0]) + SQR(E->Get(x-1,y+1,z,0,a)[0])) * NormFactor);
							TYYmod+=((1.0 / (4.0 * SQR(U->a[0]) * SQR(U->a[1]))) * (SQR(E->Get(x,y+1,z,2,a)[0]) + SQR(E->Get(x,y+1,z-1,2,a)[0])) * NormFactor);
							TYYmod+=-((1.0 / (SQR(U->a[0]) * SQR(U->a[2]))) * 0.5 * (E1Up[a] * E->Get(x,y+1,z,1,a)[0]) * NormFactor);

							//TZZmod Electric part//
							TZZmod+=((1.0 / (4.0 * SQR(U->a[1]) * SQR(U->a[2])) * (SQR(E->Get(x,y,z+1,0,a)[0]) + SQR(E->Get(x-1,y,z+1,0,a)[0]))) * NormFactor);
							TZZmod+=((1.0 / (4.0 * SQR(U->a[0]) * SQR(U->a[2])) * (SQR(E->Get(x,y,z+1,1,a)[0]) + SQR(E->Get(x,y-1,z+1,1,a)[0]))) * NormFactor);
							TZZmod+=-((1.0 / (SQR(U->a[0]) * SQR(U->a[1]))) * 0.5 * (E2Up[a] * E->Get(x,y,z+1,2,a)[0]) * NormFactor);

							//T0Xmod
							T0Xmod+=- (OJcoef[0] * E->Get(x,y,z,1,a)[0] * ((-1)*ReTrITaUxy[a]) * NormFactor );
							T0Xmod+=- (OJcoef[0] * E->Get(x+1,y,z,1,a)[0] * ((-1)*ReTrITaDx_Uxy_Ux[a]) * NormFactor );
							T0Xmod+=- (OJcoef[0] * E->Get(x,y-1,z,1,a)[0] * ((-1)*ReTrITaUymy_UxMy_Dymy[a]) * NormFactor );
							T0Xmod+=- (OJcoef[0] * E->Get(x+1,y-1,z,1,a)[0] * ((-1)*ReTrITaDxmyUymy_UxMy_DymyUxmy[a]) * NormFactor );

							T0Xmod+=- (OJcoef[0] * E->Get(x,y,z,2,a)[0] * ReTrITaUzx[a] * NormFactor );
							T0Xmod+=- (OJcoef[0] * E->Get(x+1,y,z,2,a)[0] * ReTrITaDx_Uzx_Ux[a] * NormFactor );
							T0Xmod+=- (OJcoef[0] * E->Get(x,y,z-1,2,a)[0] * ReTrITaUzmz_UMzx_Dzmz[a] * NormFactor );
							T0Xmod+=- (OJcoef[0] * E->Get(x+1,y,z-1,2,a)[0] * ReTrITaDxmzUzmz_UMzx_DzmzUxmz[a] * NormFactor );

							//T0Ymod
							T0Ymod+=- (OJcoef[1] * E->Get(x,y,z,0,a)[0] * ReTrITaUxy[a] * NormFactor );
							T0Ymod+=- (OJcoef[1] * E->Get(x,y+1,z,0,a)[0] * ReTrITaDy_Uxy_Uy[a] * NormFactor );
							T0Ymod+=- (OJcoef[1] * E->Get(x-1,y,z,0,a)[0] * ReTrITaUxmx_UMxy_Dxmx[a] * NormFactor );
							T0Ymod+=- (OJcoef[1] * E->Get(x-1,y+1,z,0,a)[0] * ReTrITaDymxUxmx_UMxy_DxmxUymx[a] * NormFactor );

							T0Ymod+=- (OJcoef[1] * E->Get(x,y,z,2,a)[0] * ((-1)*ReTrITaUyz[a]) * NormFactor );
							T0Ymod+=- (OJcoef[1] * E->Get(x,y+1,z,2,a)[0] * ((-1)*ReTrITaDy_Uyz_Uy[a]) * NormFactor );
							T0Ymod+=- (OJcoef[1] * E->Get(x,y,z-1,2,a)[0] * ((-1)*ReTrITaUzmz_UyMz_Dzmz[a]) * NormFactor );
							T0Ymod+=- (OJcoef[1] * E->Get(x,y+1,z-1,2,a)[0] * ((-1)*ReTrITaDymzUzmz_UyMz_DzmzUymz[a]) * NormFactor );

							//T0Zmod
							T0Zmod+=- (OJcoef[2] * E->Get(x,y,z,0,a)[0] * ((-1)*ReTrITaUzx[a]) * NormFactor );
							T0Zmod+=- (OJcoef[2] * E->Get(x,y,z+1,0,a)[0] * ((-1)*ReTrITaDz_Uzx_Uz[a]) * NormFactor );
							T0Zmod+=- (OJcoef[2] * E->Get(x-1,y,z,0,a)[0] *((-1)* ReTrITaUxmx_UzMx_Dxmx[a]) * NormFactor );
							T0Zmod+=- (OJcoef[2] * E->Get(x-1,y,z+1,0,a)[0] *((-1)* ReTrITaDzmxUxmx_UzMx_DxmxUzmx[a]) * NormFactor );

							T0Zmod+=- (OJcoef[2] * E->Get(x,y,z,1,a)[0] * ReTrITaUyz[a] * NormFactor );
							T0Zmod+=- (OJcoef[2] * E->Get(x,y,z+1,1,a)[0] * ReTrITaDz_Uyz_Uz[a] * NormFactor );
							T0Zmod+=- (OJcoef[2] * E->Get(x,y-1,z,1,a)[0] * ReTrITaUymy_UMyz_Dymy[a] * NormFactor );
							T0Zmod+=- (OJcoef[2] * E->Get(x,y-1,z+1,1,a)[0] * ReTrITaDzmyUymy_UMyz_DymyUzmy[a] * NormFactor );
							
						}
							

					}	
				}
			}
			
				//WRITE OUTPUT TO FILE
				OutStream << z <<"  "<< T00 <<"  "<< TXX <<"  "<< TYY <<"  "<< TZZ; //1-5
				OutStream <<"  "<< T0X <<"  "<< T0Y <<"  "<< T0Z; //6-8
				OutStream <<"  "<< T0Xmod <<"  "<< T0Ymod <<"  "<< T0Zmod; //9-11
				OutStream <<"  "<< TXXmod <<"  "<< TYYmod <<"  "<< TZZmod; //12-14
				//OutStream <<"  "<< TXY <<"  "<< TYZ <<"  "<< TZX; //15-17
			
				OutStream<<std::endl;
			
        }

        //CLOSE OUTPUT STREAM
        OutStream.close();
			
		}
		
		
		//CREATE ENERGY DENSITY MAP
		void CreateMap(std::string fname, GaugeLinks *U,GaugeLinks *UOld,ElectricFields *E){
			
			CreateMap(fname, 0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,UOld,E);
			
		}
		
		void CreateMap(std::string fname){
			
			CreateMap(fname, SUNcGaugeLinks::U,SUNcGaugeLinks::UOld,SUNcElectricFields::E);
			
		}
		
    }
}
#endif
