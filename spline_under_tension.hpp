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




#ifndef __SPLINE_U_TENSION_HPP
#define __SPLINE_U_TENSION_HPP

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


	class spline_under_tension{

		private:
			float sigmap;
			/// x values
			std::vector<float> xValues;
			/// y values
			std::vector<float> yValues; 
			/// second derivative
			std::vector<float> secondDerivative;
			int size;

			/// 
			int getInterval(float value, std::vector<float>& array, int size);
			/// 
			int binarySearch(std::vector<float>& input, float value, int low, int high, int lastMid);
			
			/// use in the spline under tension calculation
			class coefficient;
		
			/// the second derivative calculation
			class acceleration;
		public:	

			spline_under_tension();

			
			spline_under_tension(float sigma, 
					std::vector<float>& xVal,
					std::vector<float>& yVal,
					float slp1,
					float slp2,
					int islpsw);

			/// get the interpolated value
			float getInterpValue(float value);

			std::vector<float> getSecondDer(){

				return secondDerivative;
			}
			
			~spline_under_tension();









	};


#endif
