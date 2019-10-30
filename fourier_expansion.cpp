#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <functional>
#include <iomanip>
#include <fstream>

double simpson(double a, double fa, double b, double fb, double fmid)
{
	return (b - a) / 6.0 * (fa + 4.0 * fmid + fb);
}

double adpSimpsonRecursive(std::function<double(double)> f, double a, double fa, double b, double fb, double mid, double fmid, double whole, double eps, int depth)
{
	double lmid((a + mid) / 2.0), flmid(f(lmid));
	double rmid((mid + b) / 2.0), frmid(f(rmid));
	double left(simpson(a, fa, mid, fmid, flmid));
	double right(simpson(mid, fmid, b, fb, frmid));
	double delta(left + right - whole);
	
	const double epsMin(1e-9);
	const int depthMin(6);
	if (depth >= depthMin && (eps < epsMin || std::abs(delta) < 15.0 * eps))
	{
		return left + right + delta / 15.0;
	}
	else
	{
		return   adpSimpsonRecursive(f, a, fa, mid, fmid, lmid, flmid, left, eps / 2.0, depth + 1)
			   + adpSimpsonRecursive(f, mid, fmid, b, fb, rmid, frmid, right, eps / 2.0, depth + 1);
	}
}

double adpSimpson(std::function<double(double)> f, double a, double b, double eps)
{
	double fa(f(a)), fb(f(b)), mid((a + b) / 2.0), fmid(f(mid));
	return adpSimpsonRecursive(f, a, fa, b, fb, mid, fmid, simpson(a, fa, b, fb, fmid), eps, 0);
}

const double pi(3.14159265);

std::pair<std::vector<double>, std::vector<double> > getFourierSeries(std::function<double(double)> f, std::pair<double, double> domain, int samples, double eps)
{
	std::pair<std::vector<double>, std::vector<double> > series;
	double P(domain.second - domain.first);
	series.first.push_back(2.0 / P * adpSimpson(f, domain.first, domain.second, eps));
	series.second.push_back(0.0);
	for (int n(1); n < samples; ++n)
	{
		series.first.push_back(2.0 / P * adpSimpson([f, n, P](double x) -> double { return f(x) * std::cos(2.0 * pi * x * n / P); }, domain.first, domain.second, eps));
		series.second.push_back(2.0 / P * adpSimpson([f, n, P](double x) -> double { return f(x) * std::sin(2.0 * pi * x * n / P); }, domain.first, domain.second, eps));
	}
	
	return series;
}

void printFourierSeries(std::ofstream &out, std::function<double(double)> f, std::pair<double, double> domain, int samples, double eps, char symbol)
{
	
	double P{domain.second - domain.first};
	auto series(getFourierSeries(f, domain, samples, eps));
	out << std::fixed << std::setprecision(-std::log10(eps));
	
	if (std::abs(series.first[0] / 2.0) > 0.5 * eps)
		out << std::showpos << series.first[0] / 2.0;
	for (int i(1); i < samples; ++i)
	{
		double coeff(2.0 * pi * i / P);
		if (std::abs(series.first[i]) > 0.5 * eps)
			out << std::showpos << series.first[i] << "*cos(" << std::noshowpos << coeff << "*" << symbol << ")";
		if (std::abs(series.second[i]) > 0.5 * eps)
			out << std::showpos << series.second[i] << "*sin(" << std::noshowpos << coeff << "*" << symbol << ")";
	}
}

double xt(double t)
{
	const double unit{1.0 / 8.0};
	if (t <= unit * 3.0)
		return unit * 0.0 + t;
	else if (t <= unit * 7.0)
		return unit * 3.0;
	else if (t <= unit * 9.0)
		return unit * 3.0 - (t - unit * 7.0);
	else if (t <= unit * 12.0)
		return unit * 1.0;
	else if (t <= unit * 15.0)
		return unit * 1.0 + (t - unit * 12.0);
	else if (t <= unit * 22.0)
		return unit * 4.0;
	else if (t <= unit * 26.0)
		return unit * 4.0 + (t - unit * 22.0);
	else if (t <= unit * 32.0)
		return unit * 8.0;
	else if (t <= unit * 33.0)
		return unit * 8.0 - (t - unit * 32.0);
	else if (t <= unit * 38.0)
		return unit * 7.0;
	else if (t <= unit * 40.0)
		return unit * 7.0 - (t - unit * 38.0);
	else if (t <= unit * 46.0)
		return unit * 5.0;
	else if (t <= unit * 49.0)
		return unit * 5.0 + (t - unit * 46.0);
	else if (t <= unit * 50.0)
		return unit * 8.0;
	else if (t <= unit * 58.0)
		return unit * 8.0 - (t - unit * 50.0);
	else if (t <= unit * 63.0)
		return unit * 0.0;
	else if (t <= unit * 65.0)
		return unit * 0.0 + (t - unit * 63.0);
	else if (t <= unit * 67.0)
		return unit * 2.0;
	else if (t <= unit * 69.0)
		return unit * 2.0 - (t - unit * 67.0);
	else
		return unit * 0.0;
}

double yt(double t)
{
	const double unit{1.0 / 8.0};
	if (t <= unit * 3.0)
		return unit * 8.0;
	else if (t <= unit * 7.0)
		return unit * 8.0 - (t - unit * 3.0);
	else if (t <= unit * 9.0)
		return unit * 4.0;
	else if (t <= unit * 12.0)
		return unit * 4.0 - (t - unit * 9.0);
	else if (t <= unit * 15.0)
		return unit * 1.0;
	else if (t <= unit * 22.0)
		return unit * 1.0 + (t - unit * 15.0);
	else if (t <= unit * 26.0)
		return unit * 8.0;
	else if (t <= unit * 32.0)
		return unit * 8.0 - (t - unit * 26.0);
	else if (t <= unit * 33.0)
		return unit * 2.0;
	else if (t <= unit * 38.0)
		return unit * 2.0 + (t - unit * 33.0);
	else if (t <= unit * 40.0)
		return unit * 7.0;
	else if (t <= unit * 46.0)
		return unit * 7.0 - (t - unit * 40.0);
	else if (t <= unit * 49.0)
		return unit * 1.0;
	else if (t <= unit * 50.0)
		return unit * 1.0 - (t - unit * 49.0);
	else if (t <= unit * 58.0)
		return unit * 0.0;
	else if (t <= unit * 63.0)
		return unit * 0.0 + (t - unit * 58.0);
	else if (t <= unit * 65.0)
		return unit * 5.0;
	else if (t <= unit * 67.0)
		return unit * 5 + (t - unit * 65.0);
	else if (t <= unit * 69.0)
		return unit * 7.0;
	else
		return unit * 7.0 + (t - unit * 69.0);
}

int main()
{
	std::ofstream out("out.txt");
	out << "x=";
	//out << "(";
	printFourierSeries(out, xt, std::make_pair(0.0, 70.0 / 8.0), 25, 0.0001, 't');
	out << "\n\ny=";
	//out << ",";
	printFourierSeries(out, yt, std::make_pair(0.0, 70.0 / 8.0), 25, 0.0001, 't');
	out << '\n';
	//out << ")";
	return 0;
}

