/**
 * @file spline_under_tension.hpp
 * @author  Bill Brouwer <whiskeyjulietb@gmail.com>
 * @version 1.0
 *
 * @section LICENSE
 * Copyright 2013 William J. Brouwer
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * @section DESCRIPTION
 *
 * A port of the spline under tension package of A. K. CLINE AND R. J. RENKA
 *
 */



#include "spline_under_tension.hpp"

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>



class spline_under_tension::coefficient{


	private:
		float c1;
		float c2;
		float c3;

	public:
		coefficient(){}

		coefficient(float del1,
				float del2,
				float sigma,
				bool calculateAll){

			if (!calculateAll){

				c1 = -1.0f / del1;
				c2 = -c1;
			} else {


				if ((int) sigma == 0){
					//no tension
					float del = del2-del1;
					c1 = -(del1+del2) / (del1*del2);
					c2 = del2 / (del1*del);
					c3 = -del1 / (del2*del);


				} else {
					//tension
					float coshm1 = coshf(sigma*del1) -1.0f;
					float coshm2 = coshf(sigma*del2) -1.0f;

					float delp = (sigma * (del2+del1) / 2.0f);
					float delm = (sigma * (del2-del1) / 2.0f);

					float sinhmp = sinhf(delp);
					float sinhmm = sinhf(delm);

					float denom = coshm1 * (del2-del1) - 2.0f * del1 * sinhmp * sinhmm;

					c1 = 2.0f * sinhmp * sinhmm / denom;
					c2 = -coshm2 / denom;
					c3 = coshm1 /denom;

				}
			}


		}

		~coefficient(){


		}

		float getCoeffOne(){ return c1;}
		float getCoeffTwo(){ return c2;}
		float getCoeffThree(){ return c3;}

};



class spline_under_tension::acceleration{

	private:
		std::vector<float> yOutput;
		float* diagTermsPrior;
		float* diagTermsNext;
		float slpp1, slppn;

	public:
		acceleration(){}
		acceleration(float sigmap,
				std::vector<float>& xValues,
				std::vector<float>& yValues,
				float slp1,
				float slpn,
				int islpsw){


			int n=xValues.size();

			if (n<=1)
				sparseDataError();
			if (xValues[n-1] <= xValues[0])
				monotonicError();

			float delx1, delx2, delxn, delxnm;

			coefficient firstConstants, lastConstants;


			switch (islpsw){

				case 0:

					slpp1=slp1;
					slppn=slpn;
					break;

				case 1:

					slpp1=slp1;
					delxn=xValues[n-1]-xValues[n-2];
					delxnm=delxn+delxn;

					if (delxn <= 0 || delxnm <= delxn)
						monotonicError();

					lastConstants = spline_under_tension::coefficient(-delxn,-delxnm,sigmap,false);
					slppn = lastConstants.getCoeffOne() * yValues[n-1] + lastConstants.getCoeffTwo() * yValues[n-2];

					break;

					//last given
				case 2:

					slppn=slpn;
					delx1=xValues[1]-xValues[0];
					delx2=delx1+delx1;

					if (delx1 <= 0 || delx2 <= delx1)
						monotonicError();


					firstConstants = spline_under_tension::coefficient(delx1,delx2,sigmap,false);
					slpp1 = firstConstants.getCoeffOne() * yValues[0] + firstConstants.getCoeffTwo() * yValues[1];

					break;

					//both calculated
				case 3:
					delx1 = xValues[1]-xValues[0];
					delx2=xValues[2]-xValues[0];

					delxn=xValues[n-1]-xValues[n-2];
					delxnm=xValues[n-1]-xValues[n-3];

					firstConstants = spline_under_tension::coefficient(delx1,delx2,sigmap,true);
					lastConstants = spline_under_tension::coefficient(-delxn,-delxnm,sigmap,true);

					slpp1 = firstConstants.getCoeffOne() * yValues[0] + firstConstants.getCoeffTwo() * yValues[1] +
						firstConstants.getCoeffThree() * yValues[2];

					slppn = lastConstants.getCoeffOne() * yValues[n-1] + lastConstants.getCoeffTwo() * yValues[n-2] +
						lastConstants.getCoeffThree() * yValues[n-3];

					break;                                                                                        


			}


			secondDerivative(xValues,yValues,sigmap);

		}

		void secondDerivative(std::vector<float>& xValues, std::vector<float>& yValues, float sigmap){

			int n = xValues.size();

			float delx1,delx2,dx1,dx2,diag;
			for (int i=0; i<n; i++) yOutput.push_back(-1);

			delx1 = xValues[1]-xValues[0];
			std::vector<float> temp(n); 

			if (delx1 == 0)
				monotonicError();

			dx1 = (yValues[1]-yValues[0]) / delx1;

			diagTermsPrior = tridiagonalTerms(sigmap,delx1);
			yOutput[0] = (dx1-slpp1)/diagTermsPrior[0];
			temp[0] = diagTermsPrior[1] / diagTermsPrior[0];


			if (n > 2){

				for (int i=1; i<n-1; i++){

					delx2 = xValues[i+1]-xValues[i];


					if (delx2 == 0)
						monotonicError();

					dx2 = (yValues[i+1]-yValues[i]) / delx2;
					diagTermsNext = tridiagonalTerms(sigmap,delx2);
					diag = diagTermsPrior[0] + diagTermsNext[0] - diagTermsPrior[1]*temp[i-1];
					yOutput[i] = (dx2-dx1-diagTermsPrior[1]*yOutput[i-1]) / diag;

					temp[i] = diagTermsNext[1] / diag;
					dx1=dx2;

					diagTermsPrior[0]=diagTermsNext[0];
					diagTermsPrior[1]=diagTermsNext[1];


				}

			}

			diag = diagTermsPrior[0] - diagTermsPrior[1]*temp[n-2];
			yOutput[n-1] = (slppn-dx1-diagTermsPrior[1]*yOutput[n-2]) / diag;


			for (int i=n-2; i>0; i--){

				yOutput[i] = yOutput[i] - temp[i]*yOutput[i+1];

			}

		}

		float* tridiagonalTerms(float sigma, float del){


			float* ret = new float[2];

			if (sigma==0){

				ret[0] = del / 3.0f;
				ret[1] = del / 6.0f;

			} else{

				double sigdel = (sigma * del);
				double denom = del / (sinhf(sigdel) * sigdel * sigdel);

				ret[0] = (float) (denom * (sigdel * (coshf(sigdel)-1) - sinhf(sigdel) + sigdel));
				ret[1] = (float) (denom * (sinhf(sigdel) - sigdel));
			}		


			return ret;

		}


		std::vector<float> getOutput(){


			return yOutput;
		}		

		void sparseDataError(){

			std::cerr << "Insufficient samples in data \n";
			exit(1);
		}

		void monotonicError(){

			std::cerr << "Data not strictly increasing \n";
			exit(1);
		}

		~acceleration(){
		}

};
float spline_under_tension::getInterpValue(float value){


	int im1 = getInterval(value, xValues, size);

	float del1 = value - xValues[im1];
	float del2 = xValues[im1+1] - value;
	float dels = xValues[im1+1] - xValues[im1];
	float sum = (yValues[im1+1]*del1  + yValues[im1]*del2 ) / dels;

	if ((int) sigmap == 0){
		return sum - del1*del2*(secondDerivative[im1+1]*(del1+dels)
			+ secondDerivative[im1]* (del2+dels)) / (6.0f * dels);

	} else {


		float delp1 = sigmap * (del1+dels) / 2.0f;
		float delp2 = sigmap * (del2+dels) / 2.0f;

		float sinhm1 = sinhf(sigmap*del1) - sigmap*del1;
		float sinhm2 = sinhf(sigmap*del2) - sigmap*del2;
		float sinhms = sinhf(sigmap*dels) - sigmap*dels;

		float sinhp1 = sinhf(sigmap*del1 / 2.0f) - sigmap*del1 / 2.0f;
		float sinhp2 = sinhf(sigmap*del2 / 2.0f) - sigmap*del2 / 2.0f;

		float coshp1 = coshf(delp1)-1;
		float coshp2 = coshf(delp2)-1;

		return sum + (secondDerivative[im1+1] * (sinhm1*del2-del1*(2.0f * (coshp1+1.0f) * sinhp2 + sigmap*coshp1*del2))
				+ secondDerivative[im1] * (sinhm2*del1 - del2*(2.0f * (coshp2+1.0f) * sinhp1 + sigmap*coshp2*del1)))
			/ (sigmap * sigmap *dels * (sinhms + sigmap*dels));

	}

}

int spline_under_tension::getInterval(float value, std::vector<float>& array, int size){

	return binarySearch(array, value, 0, array.size()-1, 0);

}

int spline_under_tension::binarySearch(std::vector<float>& input, float value, int low, int high, int lastMid){

	int mid;

	if ((high < low) && (lastMid != (input.size() / 2)))
		return (lastMid % 2 == 0) ? lastMid : lastMid +1;
	else if ((high < low) && (lastMid == (input.size() / 2)))
		return -1;

	mid = low + (high-low) / 2;

	if (input[mid] && input[mid+1] > value)
		return binarySearch(input, value, low, mid-1, mid);
	else if (input[mid] && input[mid+1] < value)
		return binarySearch(input, value, mid+1, high, mid);
	else
		return mid;



};


spline_under_tension::spline_under_tension(){};

spline_under_tension::spline_under_tension(float sigma, 
		std::vector<float>& xVal,
		std::vector<float>& yVal,
		float slp1,
		float slp2,
		int islpsw){

	int n = xVal.size();

	sigmap = fabs(sigma) * (float) (n-1) / (xVal[n-1]-xVal[0]);
	xValues = xVal;
	yValues = yVal;
	size = n;

	secondDerivative = acceleration(sigmap,xValues,yValues,slp1,slp2,islpsw).getOutput();
}


spline_under_tension::~spline_under_tension(){};

