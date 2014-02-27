#include <iostream>
#include <iterator>
#include "spline_under_tension.hpp"
using namespace std;


int main(){


	vector<float> xvals,yvals;

	for (int i=0; i<2048; i++){
		xvals.push_back(i*0.05f);
		yvals.push_back(sin(xvals[i]));
	}


	spline_under_tension foo(0.1, xvals, yvals, -1, -1, 3);

	//vector<float> out=foo.getSecondDer();

	//for (vector<float>::iterator it = out.begin(); it != out.end(); it++)
	//	cout << *it << endl;

	//exit(0);

	cout << "#yvals; interp" << endl;
	cout << foo.getInterpValue(xvals[0]+0.03333);
	for (int i=1; i<1024; i++){
		cout << "," << foo.getInterpValue(xvals[i]+0.3333);

	}

	cout << endl;
	cout << "#actual yvals" << endl;
	cout << sin(xvals[0]+0.03333);
	for (int i=1; i<1024; i++){
		cout << "," << sin(xvals[i]+0.3333);
	}
	cout << endl;
	cout << "#xvals" << endl;
	cout << xvals[0]+0.03333;
	for (int i=1; i<1024; i++){
		cout << ","<< xvals[i]+0.3333;

	}
	cout << endl;

	return 0;
};
