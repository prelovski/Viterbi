
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <random>
#include <sstream>

using namespace std;


using bit = unsigned char;
#define ZERO 0x0
#define ONE 0x1
using poly = unsigned long long;

int popcnt(poly p)
{
	int cnt = 0;
	while (p!=0)
	{
		if ((p & 1) == 0)
			cnt++;
		p >>= 1;
	}
	return cnt;
}

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "Russian");

	const double SNR = atoi(argv[1]);
	const int k = atoi(argv[2]);
	const int np = atoi(argv[3]);
	vector<poly> p(np);
	
	stringstream ss;
	for (int i = 0; i < np; i++)
		ss << argv[4 + i] << ' ';
	for (int i = 0; i < np; i++)
		ss >> oct >> p[i];

	//длина кодового ограничени€ - максимальна€ степень полинома
	int W = 0;
	for (int i = 0; i < np; i++)
		W = max(W, popcnt(p[i]));


	vector<bit> u(k+W-1,ZERO);
	vector<bit> c((k + W - 1)*np, ZERO);
	
	mt19937 genRandom;
	uniform_int_distribution<> dis(0, 1);
	for (int i = W - 1; i < u.size(); i++)
		u[i] = dis(genRandom);

	//вывод вектора u на экран
	cout << "u = ";
	for (int i = 0; i < u.size(); ++i)
		cout << int(u[i]);
	cout << endl;

	for (int i = 0; i < p.size(); i++){
		poly p_i = p[i];
		int shift = 0;
		while (p_i != 0) {
			if ((p_i & 1) != 0) {
				for (int j = shift; j < u.size(); j++){
					c[i*u.size()+j - shift] ^= u[j];
				}
			}
			p_i >>= 1, shift++;
		}
	}

	//вывод вектора u на экран
	cout << "с = ";
	for (int i = 0; i < c.size(); ++i)
		cout << int(c[i]);
	cout << endl;

	//определ€ем параметры дл€ гауссовского канала (на основании SNR - отношение сигнал\шум и R - скорости кода)
	double R = 1.0 / np;
	double sigma = sqrt(pow(10, -SNR / 10) / (2 * R));
	normal_distribution<> ndis(0, sigma);

	cout << "sigma= " << sigma<<endl;

	vector<double> noisy((k + W - 1)*np, ZERO);
	vector<double> LLR((k + W - 1)*np, ZERO);
	for (int i = 0; i < noisy.size(); i++) {
		noisy[i] = (c[i]==ZERO?1:-1) + ndis(genRandom);
		LLR[i] = 2 * noisy[i] / (sigma*sigma);
	}



	//вывод вектора noisy на экран
	cout << "noisy = ";
	for (int i = 0; i < noisy.size(); ++i)
		cout << noisy[i] << " ";
	cout << endl;


	//вывод вектора LLR на экран
	cout << "LLR = ";
	for (int i = 0; i < LLR.size(); ++i)
		cout << LLR[i] << " ";
	cout << endl;

	//вычисл€ем Ћќѕы: l[i] = ln( P(u[i] = 0 | y[i]) / P(u[i] = 1 | y[i]) )
	/*for (int i = 0; i < noisy.size(); ++i) {
		noisy[i] = 2*
	}*/

	//system("PAUSE");
	return 0;
}