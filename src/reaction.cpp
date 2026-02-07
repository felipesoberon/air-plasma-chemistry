#include "reaction.h"

// Static BOLSIG+ table member definitions
bool reaction::bolsigLoaded = false;
int reaction::bolsigNPoints = 0;
double reaction::bolsigTeEv[reaction::bolsigMaxPoints] = {};
double reaction::bolsigRate[reaction::bolsigMaxPoints][reaction::bolsigMaxReactions] = {};
int reaction::bolsigNCols = 0;
int reaction::bolsigPengCol[674] = {};  // initialized in loadBolsigTable


reaction::reaction(void)
{
	int i;
	for (i=0; i<=maxNoReactionSpecies; i++)
    {
		reactantSpeciesList[i] = 0;
		productSpeciesList[i] = 0;
	}
}



void reaction::loadBolsigTable(const char* filename)
{
	// Initialize all mapping entries to -1 (not present)
	for (int i = 0; i < 674; i++) bolsigPengCol[i] = -1;

	ifstream file(filename);
	if (!file.is_open())
	{
		cout << "WARNING: Could not open BOLSIG+ rates file: " << filename << endl;
		bolsigLoaded = false;
		return;
	}

	// Parse header row
	string headerLine;
	getline(file, headerLine);
	stringstream hss(headerLine);
	string token;
	int colIndex = 0;
	bolsigNCols = 0;

	while (getline(hss, token, ','))
	{
		// Remove whitespace/carriage returns
		while (!token.empty() && (token.back() == '\r' || token.back() == ' '))
			token.pop_back();

		// Check if token starts with 'R' followed by digits
		if (token.length() > 1 && token[0] == 'R')
		{
			int pengNum = 0;
			bool valid = true;
			// Handle tokens like "R653_C18" by only parsing digits after 'R'
			size_t k = 1;
			while (k < token.length() && token[k] >= '0' && token[k] <= '9') k++;
			if (k > 1)
			{
				pengNum = atoi(token.substr(1, k-1).c_str());
				if (pengNum >= 625 && pengNum <= 673)
				{
					// Only store the first occurrence (skip duplicates like R653_C18)
					if (bolsigPengCol[pengNum] == -1)
					{
						bolsigPengCol[pengNum] = bolsigNCols;
						bolsigNCols++;
					}
					else
					{
						// This is a duplicate column (e.g. R653_C18 after R653)
						// Store it as additional column but don't overwrite mapping
						bolsigNCols++;
					}
				}
			}
		}
		colIndex++;
	}

	// Re-parse header to build a column-to-bolsigRate-index mapping
	// We need to know which data column maps to which bolsigRate column
	int totalDataCols = colIndex;
	int* dataColToBolsigCol = new int[totalDataCols];
	for (int i = 0; i < totalDataCols; i++) dataColToBolsigCol[i] = -1;

	// Reset and re-parse
	for (int i = 0; i < 674; i++) bolsigPengCol[i] = -1;
	bolsigNCols = 0;

	stringstream hss2(headerLine);
	colIndex = 0;
	while (getline(hss2, token, ','))
	{
		while (!token.empty() && (token.back() == '\r' || token.back() == ' '))
			token.pop_back();

		if (token.length() > 1 && token[0] == 'R')
		{
			size_t k = 1;
			while (k < token.length() && token[k] >= '0' && token[k] <= '9') k++;
			if (k > 1)
			{
				int pengNum = atoi(token.substr(1, k-1).c_str());
				if (pengNum >= 625 && pengNum <= 673)
				{
					if (bolsigPengCol[pengNum] == -1)
					{
						bolsigPengCol[pengNum] = bolsigNCols;
					}
					dataColToBolsigCol[colIndex] = bolsigNCols;
					bolsigNCols++;
				}
			}
		}
		colIndex++;
	}

	// Read data rows
	bolsigNPoints = 0;
	string dataLine;
	while (getline(file, dataLine))
	{
		if (dataLine.empty()) continue;
		if (bolsigNPoints >= bolsigMaxPoints) break;

		stringstream dss(dataLine);
		string val;
		int col = 0;
		while (getline(dss, val, ','))
		{
			while (!val.empty() && (val.back() == '\r' || val.back() == ' '))
				val.pop_back();

			if (col == 1)  // Te_eV is the second column
			{
				bolsigTeEv[bolsigNPoints] = atof(val.c_str());
			}
			else if (col >= 2 && col < totalDataCols)
			{
				int bCol = dataColToBolsigCol[col];
				if (bCol >= 0 && bCol < bolsigMaxReactions)
				{
					bolsigRate[bolsigNPoints][bCol] = atof(val.c_str());
				}
			}
			col++;
		}
		bolsigNPoints++;
	}

	delete[] dataColToBolsigCol;
	file.close();
	bolsigLoaded = true;
	cout << "NOTICE: Loaded BOLSIG+ rates table: " << bolsigNPoints << " points, " << bolsigNCols << " reaction columns." << endl;
}



double reaction::interpolateBolsigRate(int pengReaction, double TeEv)
{
	if (!bolsigLoaded) return 0.0;
	if (pengReaction < 625 || pengReaction > 673) return 0.0;

	int col = bolsigPengCol[pengReaction];
	if (col < 0) return 0.0;

	// Clamp to range
	if (TeEv <= bolsigTeEv[0])
		return bolsigRate[0][col];
	if (TeEv >= bolsigTeEv[bolsigNPoints - 1])
		return bolsigRate[bolsigNPoints - 1][col];

	// Find bracketing interval
	for (int i = 0; i < bolsigNPoints - 1; i++)
	{
		if (TeEv >= bolsigTeEv[i] && TeEv <= bolsigTeEv[i + 1])
		{
			double t = (TeEv - bolsigTeEv[i]) / (bolsigTeEv[i + 1] - bolsigTeEv[i]);
			return bolsigRate[i][col] + t * (bolsigRate[i + 1][col] - bolsigRate[i][col]);
		}
	}

	return bolsigRate[bolsigNPoints - 1][col];
}



double reaction::reactionRateFunction(int j, double Tgas, double Telectron)
{
	double result = 0;
	double const kb = 1.38064852E-23;  //Boltzmann constant 
	double const qe = 1.6021766208E-19; //Electron Charge
	double eps; 
	double Te = Telectron;
	double Tg = Tgas;
	eps = 3 * kb * Te / ( 2 * qe );  //electron mean energy
	
	switch(j)
    {
		case 1 :
		if (Te!=0) result = 3.99E-17 * pow(eps,2.24) * exp(-9.1/eps);
		break;
		
		case 2 :
		if (Te!=0) result = 3.34E-16 * pow(eps,-0.06) * exp(-8.5/eps);
		break;
		
		case 3 :
		if (Te!=0) result = 8.44E-15 * pow(eps,-0.33) * exp(-9.15/eps);
		break;
		
		case 4 :
		if (Te!=0) result = 1E-16 * pow(eps,1.90) * exp(-14.6/eps);
		break;
		
		case 5 :
		if (Te!=0) result = 5.06E-15 * exp(-10.8/( pow(eps,3.95) ));
		break;
		
		case 6 :
		if (Te!=0) result = 1.45E-17 * pow(eps,2.58) * exp(-8.54/eps);
		break;
		
		case 7 :
		if (Te!=0) result = 2.03E-14 * pow(eps,-0.10) * exp(-8.47/eps);
		break;
		
		case 8 :
		if (Te!=0) result = 1.82E-14 * pow(eps,-0.13) * exp(-10.7/eps);
		break;
		
		case 9 :
		if (Te!=0) result = 1.04E-15 * exp(-2.59/eps);
		break;
		
		case 10 :
		if (Te!=0) result = 9.54E-12 * pow(eps,-1.05) * exp(-55.6/eps);
		break;
		
		case 11 :
		if (Te!=0) result = 1.78E-12 * pow(eps,-0.614) * exp(-11.5/eps);
		break;
		
		case 12 :
		if (Te!=0) result = 7.46E-15 * exp(-5.58/pow(eps,1.47));
		break;
		
		case 13 :
		if (Te!=0) result = 4.75E-15 * pow(eps,0.61) * exp(-22.1/eps);
		break;
		
		case 14 :
		if (Te!=0) result = 9.65E-18 * pow(eps,2.53) * exp(-8.99/eps);
		break;
		
		case 15 :
		if (Te!=0) result = 9.89E-12 * pow(eps,-1.64) * exp(-67.6/eps);
		break;
		
		case 16 :
		if (Te!=0) result = 7.45E-15 * pow(eps,0.34) * exp(-54.2/eps);
		break;
		
		case 17 :
		if (Te!=0) result = 7.4E-16 * pow(eps,0.45) * exp(-55.5/eps);
		break;
		
		case 18 :
		if (Te!=0) result = 8.49E-15 * pow(eps,-1.23) * exp(-74.0/eps);
		break;
		
		case 19 :
		if (Te!=0) result = 5.15E-15 * pow(eps,0.62) * exp(-10.9/eps);
		break;
		
		case 20 :
		if (Te!=0) result = 5.19E-18 * pow(eps,1.2) * exp(-13.8/eps);
		break;
		
		case 21 :
		if (Te!=0) result = 3.29E-15 * pow(eps,0.578) * exp(-7.56/eps);
		break;
		
		case 22 :
		if (Te!=0) result = 4E-17 * pow(eps,2.13) * exp(-14.9/eps);
		break;
		
		case 23 :
		if (Te!=0) result = 2.43E-17 * pow(eps,2.77) * exp(-5.62/eps);
		break;
		
		case 24 :
		if (Te!=0) result = 3.12E-35 / pow(Te,1.5);
		break;
		
		case 25 :
		if (Te!=0) result = 1.66E-12 / pow(Te,0.7);
		break;
		
		case 26 :
		if (Te!=0) result = 1.5E-12 / pow(Te,0.7);
		break;
		
		case 27 :
		if (Te!=0) result = 3.12E-35 / pow(Te,1.5);
		break;
		
		case 28 :
		if (Te!=0) result = 3.46E-12 / pow(Te,0.5);
		break;
		
		case 29 :
		if (Te!=0) result = 4.73E-11 / pow(Te,0.53);
		break;
		
		case 30 :
		if (Te!=0) result = 3.12E-35 / pow(Te,1.5);
		break;
		
		case 31 :
		if (Te!=0) result = 1.68E-11 / pow(Te,0.7);
		break;
		
		case 32 :
		if (Te!=0) result = 1.24E-11 / pow(Te,0.7);
		break;
		
		case 33 :
		if (Te!=0) result = 3.12E-35 / pow(Te,1.5);
		break;
		
		case 34 :
		if (Te!=0) result = 2.42E-11 / pow(Te,0.5);
		break;
		
		case 35 :
		if (Te!=0) result = 3.46E-12 / pow(Te,0.5);
		break;
		
		case 36 :
		if (Te!=0) result = 1.07E-11 / pow(Te,0.85);
		break;
		
		case 37 :
		if (Te!=0) result = 4.28E-11 / pow(Te,0.85);
		break;
		
		case 38 :
		if (Te!=0) result = 3.12E-35 / pow(Te,1.5);
		break;
		
		case 39 :
		if (Te!=0) result = 3.46E-12 / pow(Te,0.5);
		break;
		
		case 40 :
		if (Te!=0) result = 1.86E-13 / pow(Te,0.43);
		break;
		
		case 41 :
		if (Te!=0) result = 5.20E-11 / pow(Te,0.5);
		break;
		
		case 42 :
		if (Te!=0) result = 1.14E-11 / pow(Te,0.97);
		break;
		
		case 43 :
		if (Te!=0) result = 2.73E-12 / pow(Te,0.5);
		break;
		
		case 44 :
		if (Te!=0) result = 1.37E-12 / pow(Te,0.5);
		break;
		
		case 45 :
		if (Te!=0) result = 1.37E-12 / pow(Te,0.5);
		break;
		
		case 46 :
		if (Te!=0) result = 5.46E-12 / pow(Te,0.5);
		break;
		
		case 47 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 48 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 49 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 50 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 51 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 52 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 53 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 54 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 55 :
		if (Te!=0) result = 1E-31 * pow(Tg/Te,4.5);
		break;
		
		case 56 :
		result = 1E-43;
		break;
		
		case 57 :
		result = 1E-43;
		break;
		
		case 58 :
		if (Te!=0) result = 1.4E-41 * (Tg/Te) * exp(-600/Tg) * \
		exp(700*(Te-Tg)/(Te*Tg));
		break;
		
		case 59 :
		if (Te!=0) result = 1.1E-43 * pow(Tg/Te,2) * exp(-70/Tg) * \
		exp(1500*(Te-Tg)/(Te*Tg));
		break;
		
		case 60 :
		if (Te!=0) result = 2.63E-16 * pow(eps,-0.495) * exp(-5.65/eps);
		break;
		
		case 61 :
		if (eps>1.1) result = 9.72E-15 * pow(eps,-1.62) * exp(-14.2/eps);
		else result = 2.78E-20;
		break;
		
		case 62 :
		result = 1E-17;
		break;
		
		case 63 :
		result = 1E-15;
		break;
		
		case 64 :
		result = 1E-43;
		break;
		
		case 65 :
		result = 2E-16;
		break;
		
		case 66 :
		result = 8E-43;
		break;
		
		case 67 :
		result = 1E-17;
		break;
		
		case 68 :
		result = 1.5E-42;
		break;
		
		case 69 :
		result = 1E-42;
		break;
		
		case 70 :
		if (Te!=0) result = 4.42E-14 * pow(eps,-2.0) * exp(-13.39/eps);
		break;
		
		case 71 :
		if (Te!=0) result = 2.97E-15 * pow(eps,-1.56) * exp(-13.67/eps);
		break;
		
		case 72 :
		if (Te!=0) result = 9.6E-16 * pow(eps,-1.70) * exp(-13.31/eps);
		break;
		
		case 73 :
		result = 5E-14;
		break;
		
		case 74 :
		result = 2.6E-16;
		break;
		
		case 75 :
		result = 1E-18;
		break;
		
		case 76 :
		result = 2.2E-15;
		break;
		
		case 77 :
		result = 1.9E-15;
		break;
		
		case 78 :
		result = 1.4E-16;
		break;
		
		case 79 :
		result = 1E-18;
		break;
		
		case 80 :
		result = 3E-16;
		break;
		
		case 81 :
		result = 3E-16;
		break;
		
		case 82 :
		result = 2.6E-16;
		break;
		
		case 83 :
		result = 7E-16;
		break;
		
		case 84 :
		result = 5E-16;
		break;
		
		case 85 :
		result = 1.9E-18 * pow(Tg/300,0.5) * exp(-4990/Tg);
		break;
		
		case 86 :
		result = 2.1E-15;
		break;
		
		case 87 :
		result = 2.5E-15;
		break;
		
		case 88 :
		result = 1.5E-16;
		break;
		
		case 89 :
		result = 2.7E-16 * pow(Tg/300,0.5) * exp(-5590/Tg);
		break;
		
		case 90 :
		result = 6E-16;
		break;
		
		case 91 :
		result = 2E-16;
		break;
		
		case 92 :
		result = 1.4E-15;
		break;
		
		case 93 :
		result = 3E-16;
		break;
		
		case 94 :
		result = 2.3E-17;
		break;
		
		case 95 :
		result = 3E-16;
		break;
		
		case 96 :
		result = 2.40E-19;
		break;
		
		case 97 :
		result = 5E-18;
		break;
		
		case 98 :
		result = 5.1E-18;
		break;
		
		case 99 :
		result = 2.3E-19;
		break;
		
		case 100 :
		result = 1E-18;
		break;
		
		case 101 :
		result = 1E-18;
		break;
		
		case 102 :
		result = 1E-18;
		break;
		
		case 103 :
		result = 1E-18;
		break;
		
		case 104 :
		result = 1.2E-15;
		break;
		
		case 105 :
		result = 1.8E-15;
		break;
		
		case 106 :
		result = 2E-16;
		break;
		
		case 107 :
		result = 1.8E-15;
		break;
		
		case 108 :
		result = 1E-41;
		break;
		
		case 109 :
		result = 4.6E-41;
		break;
		
		case 110 :
		result = 5.5E-16;
		break;
		
		case 111 :
		result = 4.72E-16;
		break;
		
		case 112 :
		result = 8.33E-17;
		break;
		
		case 113 :
		result = 1E-18;
		break;
		
		case 114 :
		result = 3E-16;
		break;
		
		case 115 :
		result = 5E-16;
		break;
		
		case 116 :
		result = 1E-18;
		break;
		
		case 117 :
		result = 1E-41;
		break;
		
		case 118 :
		result = 2.7E-16;
		break;
		
		case 119 :
		result = 2.8E-17;
		break;
		
		case 120 :
		result = 3E-16;
		break;
		
		case 121 :
		result = 5E-16;
		break;
		
		case 122 :
		result = 3.4E-16;
		break;
		
		case 123 :
		result = 3.4E-16;
		break;
		
		case 124 :
		result = 1.19E-15;
		break;
		
		case 125 :
		result = 2.1E-16;
		break;
		
		case 126 :
		result = 1E-18;
		break;
		
		case 127 :
		result = 1E-41 * (300/Tg);
		break;
		
		case 128 :
		result = 1E-41 * (300/Tg);
		break;
		
		case 129 :
		result = 3E-16;
		break;
		
		case 130 :
		result = 6E-16;
		break;
		
		case 131 :
		result = 4E-16;
		break;
		
		case 132 :
		result = 3.9E-16;
		break;
		
		case 133 :
		result = 5E-17;
		break;
		
		case 134 :
		result = 3E-16;
		break;
		
		case 135 :
		result = 1.4E-16;
		break;
		
		case 136 :
		result = 1.8E-16 * (300/Tg);
		break;
		
		case 137 :
		result = 1E-17 * pow(300/Tg,0.5);
		break;
		
		case 138 :
		result = 5E-17;
		break;
		
		case 139 :
		result = 1E-16;
		break;
		
		case 140 :
		result = 2.3E-15;
		break;
		
		case 141 :
		result = 6.6E-17;
		break;
		
		case 142 :
		result = 2.3E-17;
		break;
		
		case 143 :
		result = 2E-17;
		break;
		
		case 144 :
		result = 4.4E-17;
		break;
		
		case 145 :
		result = 7E-17;
		break;
		
		case 146 :
		result = 7E-17;
		break;
		
		case 147 :
		result = 5E-17;
		break;
		
		case 148 :
		result = 7E-17;
		break;
		
		case 149 :
		result = 7E-17;
		break;
		
		case 150 :
		result = 2.1E-16 * exp(Tg/121);
		break;
		
		case 151 :
		result = 3E-16;
		break;
		
		case 152 :
		result = 1E-17;
		break;
		
		case 153 :
		result = 3.9E-16;
		break;
		
		case 154 :
		result = 2.5E-16;
		break;
		
		case 155 :
		result = 5E-17;
		break;
		
		case 156 :
		result = 2.5E-16;
		break;
		
		case 157 :
		result = 2.4E-16;
		break;
		
		case 158 :
		result = 3E-16 * exp(-1800/Tg);
		break;
		
		case 159 :
		result = 3E-15;
		break;
		
		case 160 :
		result = 1E-41;
		break;
		
		case 161 :
		result = 6E-41 *pow(300/Tg,2);
		break;
		
		case 162 :
		result = 1E-41;
		break;
		
		case 163 :
		result = 2.1E-17 * pow(300/Tg,0.5);
		break;
		
		case 164 :
		result = 1E-16;
		break;
		
		case 165 :
		result = 1.3E-16;
		break;
		
		case 166 :
		result = 1E-18;
		break;
		
		case 167 :
		result = 3E-18;
		break;
		
		case 168 :
		result = 6.3E-16;
		break;
		
		case 169 :
		result = 2.3E-16;
		break;
		
		case 170 :
		result = 2E-17;
		break;
		
		case 171 :
		result = 5E-16;
		break;
		
		case 172 :
		result = 1.6E-15;
		break;
		
		case 173 :
		result = 6.8E-16;
		break;
		
		case 174 :
		result = 1.7E-15;
		break;
		
		case 175 :
		result = 3.3E-16;
		break;
		
		case 176 :
		result = 3.6E-16;
		break;
		
		case 177 :
		result = 3.2E-15;
		break;
		
		case 178 :
		result = 5.5E-43 * pow(300/Tg,2.7);
		break;
		
		case 179 :
		result = 1.5E-16;
		break;
		
		case 180 :
		result = 1E-23;
		break;
		
		case 181 :
		result = 8.8E-16;
		break;
		
		case 182 :
		result = 4.6E-16;
		break;
		
		case 183 :
		result = 6.6E-16;
		break;
		
		case 184 :
		result = 1E-17;
		break;
		
		case 185 :
		result = 3E-16;
		break;
		
		case 186 :
		result = 3.3E-12 * pow(300/Tg,4) * exp(-5030/Tg);
		break;
		
		case 187 :
		result = 6.80E-16;
		break;
		
		case 188 :
		result = 3E-16;
		break;
		
		case 189 :
		result = 1.1E-42 * (300/Tg);
		break;
		
		case 190 :
		result = 1E-16;
		break;
		
		case 191 :
		result = 8E-16;
		break;
		
		case 192 :
		result = 2E-16;
		break;
		
		case 193 :
		result = 2E-18;
		break;
		
		case 194 :
		result = 1E-41;
		break;
		
		case 195 :
		result = 1.2E-15;
		break;
		
		case 196 :
		result = 3E-16;
		break;
		
		case 197 :
		result = 3.3E-17;
		break;
		
		case 198 :
		result = 1.4E-15;
		break;
		
		case 199 :
		result = 3.3E-16;
		break;
		
		case 200 :
		result = 3.5E-43 * (300/Tg);
		break;
		
		case 201 :
		result = 3.5E-16;
		break;
		
		case 202 :
		result = 1E-17;
		break;
		
		case 203 :
		result = 7E-16;
		break;
		
		case 204 :
		result = 5E-16;
		break;
		
		case 205 :
		result = 2.8E-16;
		break;
		
		case 206 :
		result = 1E-17;
		break;
		
		case 207 :
		result = 1E-17;
		break;
		
		case 208 :
		result = 1E-17;
		break;
		
		case 209 :
		result = 2E-17;
		break;
		
		case 210 :
		result = 7E-17;
		break;
		
		case 211 :
		result = 5E-16;
		break;
		
		case 212 :
		result = 8.4E-16;
		break;
		
		case 213 :
		result = 2.5E-16;
		break;
		
		case 214 :
		result = 3E-16;
		break;
		
		case 215 :
		result = 4E-16;
		break;
		
		case 216 :
		result = 1E-16 * exp(-1044/Tg);
		break;
		
		case 217 :
		result = 2.3E-16;
		break;
		
		case 218 :
		result = 1.2E-17;
		break;
		
		case 219 :
		result = 4.29E-16;
		break;
		
		case 220 :
		result = 2.21E-16;
		break;
		
		case 221 :
		result = 4.59E-17;
		break;
		
		case 222 :
		result = 2.24E-16;
		break;
		
		case 223 :
		result = 5.9E-16;
		break;
		
		case 224 :
		result = 1E-21;
		break;
		
		case 225 :
		result = 1E-41 * (300/Tg);
		break;
		
		case 226 :
		result = 2.8E-20;
		break;
		
		case 227 :
		result = 3E-16;
		break;
		
		case 228 :
		result = 3E-16;
		break;
		
		case 229 :
		result = 5E-16;
		break;
		
		case 230 :
		result = 3E-16;
		break;
		
		case 231 :
		result = 3E-16;
		break;
		
		case 232 :
		result = 2.75E-16;
		break;
		
		case 233 :
		result = 7E-16;
		break;
		
		case 234 :
		result = 2.75E-16;
		break;
		
		case 235 :
		result = 4E-18;
		break;
		
		case 236 :
		result = 5E-16;
		break;
		
		case 237 :
		result = 1.8E-17;
		break;
		
		case 238 :
		result = 4E-16;
		break;
		
		case 239 :
		result = 5E-19;
		break;
		
		case 240 :
		result = 1.6E-15;
		break;
		
		case 241 :
		result = 3E-21;
		break;
		
		case 242 :
		result = 3.8E-16;
		break;
		
		case 243 :
		result = 1.17E-15;
		break;
		
		case 244 :
		result = 1.9E-15;
		break;
		
		case 245 :
		result = 3.1E-41;
		break;
		
		case 246 :
		result = 8.2E-15;
		break;
		
		case 247 :
		result = 1.1E-15;
		break;
		
		case 248 :
		result = 2.9E-15;
		break;
		
		case 249 :
		result = 3.8E-15;
		break;
		
		case 250 :
		result = 7.83E-16;
		break;
		
		case 251 :
		result = 6.4E-16;
		break;
		
		case 252 :
		result = 2E-15;
		break;
		
		case 253 :
		result = 3.43E-15;
		break;
		
		case 254 :
		result = 3.86E-15;
		break;
		
		case 255 :
		result = 8E-16;
		break;
		
		case 256 :
		result = 3E-15;
		break;
		
		case 257 :
		result = 7E-16;
		break;
		
		case 258 :
		result = 5.9E-16;
		break;
		
		case 259 :
		result = 5.2E-16;
		break;
		
		case 260 :
		result = 1.3E-15;
		break;
		
		case 261 :
		result = 2.13E-16;
		break;
		
		case 262 :
		result = 9.7E-16;
		break;
		
		case 263 :
		result = 7E-16;
		break;
		
		case 264 :
		result = 1.59E-15;
		break;
		
		case 265 :
		result = 1.3E-15;
		break;
		
		case 266 :
		result = 9E-16;
		break;
		
		case 267 :
		result = 1.9E-15;
		break;
		
		case 268 :
		result = 1.9E-16;
		break;
		
		case 269 :
		result = 5.5E-17;
		break;
		
		case 270 :
		result = 4.3E-16;
		break;
		
		case 271 :
		result = 4.6E-16;
		break;
		
		case 272 :
		result = 1.2E-15;
		break;
		
		case 273 :
		result = 7.6E-16;
		break;
		
		case 274 :
		result = 1.7E-15;
		break;
		
		case 275 :
		result = 1.5E-18;
		break;
		
		case 276 :
		result = 5.5E-16;
		break;
		
		case 277 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 278 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 279 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 280 :
		result = 1E-13;
		break;
		
		case 281 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 282 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 283 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 284 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 285 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 286 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 287 :
		result = 1E-13;
		break;
		
		case 288 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 289 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 290 :
		result = 1E-13;
		break;
		
		case 291 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 292 :
		result = 1E-13;
		break;
		
		case 293 :
		result = 1E-13;
		break;
		
		case 294 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 295 :
		result = 1E-13;
		break;
		
		case 296 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 297 :
		result = 1E-13;
		break;
		
		case 298 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 299 :
		result = 1E-13;
		break;
		
		case 300 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 301 :
		result = 1E-13;
		break;
		
		case 302 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 303 :
		result = 1E-13;
		break;
		
		case 304 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 305 :
		result = 1E-13;
		break;
		
		case 306 :
		result = 1E-13;
		break;
		
		case 307 :
		result = 1E-13;
		break;
		
		case 308 :
		result = 1E-13;
		break;
		
		case 309 :
		result = 1E-13;
		break;
		
		case 310 :
		result = 1E-13;
		break;
		
		case 311 :
		result = 1E-13;
		break;
		
		case 312 :
		result = 1E-13;
		break;
		
		case 313 :
		result = 1E-13;
		break;
		
		case 314 :
		result = 1E-13;
		break;
		
		case 315 :
		result = 1E-13;
		break;
		
		case 316 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 317 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 318 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 319 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 320 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 321 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 322 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 323 :
		result = 1E-13;
		break;
		
		case 324 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 325 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 326 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 327 :
		result = 1E-13;
		break;
		
		case 328 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 329 :
		result = 1E-13;
		break;
		
		case 330 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 331 :
		result = 1E-13;
		break;
		
		case 332 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 333 :
		result = 1E-13;
		break;
		
		case 334 :
		result = 1E-13;
		break;
		
		case 335 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 336 :
		result = 1E-13;
		break;
		
		case 337 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 338 :
		result = 1E-13;
		break;
		
		case 339 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 340 :
		result = 1E-13;
		break;
		
		case 341 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 342 :
		result = 1E-13;
		break;
		
		case 343 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 344 :
		result = 1E-13;
		break;
		
		case 345 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 346 :
		result = 1E-13;
		break;
		
		case 347 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 348 :
		result = 1E-13;
		break;
		
		case 349 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 350 :
		result = 1E-13;
		break;
		
		case 351 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 352 :
		result = 1E-13;
		break;
		
		case 353 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 354 :
		result = 1E-13;
		break;
		
		case 355 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 356 :
		result = 1E-13;
		break;
		
		case 357 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 358 :
		result = 1E-13;
		break;
		
		case 359 :
		result = 1E-13;
		break;
		
		case 360 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 361 :
		result = 1E-13;
		break;
		
		case 362 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 363 :
		result = 1E-13;
		break;
		
		case 364 :
		result = 1E-13;
		break;
		
		case 365 :
		result = 1E-13;
		break;
		
		case 366 :
		result = 1E-13;
		break;
		
		case 367 :
		result = 1E-13;
		break;
		
		case 368 :
		result = 1E-13;
		break;
		
		case 369 :
		result = 1E-13;
		break;
		
		case 370 :
		result = 1E-13;
		break;
		
		case 371 :
		result = 1E-13;
		break;
		
		case 372 :
		result = 1E-13;
		break;
		
		case 373 :
		result = 1E-13;
		break;
		
		case 374 :
		result = 1E-13;
		break;
		
		case 375 :
		result = 1E-13;
		break;
		
		case 376 :
		result = 1E-13;
		break;
		
		case 377 :
		result = 1E-13;
		break;
		
		case 378 :
		result = 1E-13;
		break;
		
		case 379 :
		result = 1E-13;
		break;
		
		case 380 :
		result = 1E-13;
		break;
		
		case 381 :
		result = 1E-13;
		break;
		
		case 382 :
		result = 1E-13;
		break;
		
		case 383 :
		result = 1E-13;
		break;
		
		case 384 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 385 :
		result = 1E-13;
		break;
		
		case 386 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 387 :
		result = 1E-13;
		break;
		
		case 388 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 389 :
		result = 1E-13;
		break;
		
		case 390 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 391 :
		result = 1E-13;
		break;
		
		case 392 :
		result = 1E-13;
		break;
		
		case 393 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 394 :
		result = 1E-13;
		break;
		
		case 395 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 396 :
		result = 1E-13;
		break;
		
		case 397 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 398 :
		result = 1E-13;
		break;
		
		case 399 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 400 :
		result = 1E-13;
		break;
		
		case 401 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 402 :
		result = 1E-13;
		break;
		
		case 403 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 404 :
		result = 1E-13;
		break;
		
		case 405 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 406 :
		result = 1E-13;
		break;
		
		case 407 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 408 :
		result = 1E-13;
		break;
		
		case 409 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 410 :
		result = 1E-13;
		break;
		
		case 411 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 412 :
		result = 1E-13;
		break;
		
		case 413 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 414 :
		result = 1E-13;
		break;
		
		case 415 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 416 :
		result = 1E-13;
		break;
		
		case 417 :
		result = 1E-13;
		break;
		
		case 418 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 419 :
		result = 1E-13;
		break;
		
		case 420 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 421 :
		result = 1E-13;
		break;
		
		case 422 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 423 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 424 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 425 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 426 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 427 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 428 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 429 :
		result = 1E-13;
		break;
		
		case 430 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 431 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 432 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 433 :
		result = 1E-13;
		break;
		
		case 434 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 435 :
		result = 1E-13;
		break;
		
		case 436 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 437 :
		result = 1E-13;
		break;
		
		case 438 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 439 :
		result = 1E-13;
		break;
		
		case 440 :
		result = 1E-13;
		break;
		
		case 441 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 442 :
		result = 1E-13;
		break;
		
		case 443 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 444 :
		result = 1E-13;
		break;
		
		case 445 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 446 :
		result = 1E-13;
		break;
		
		case 447 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 448 :
		result = 1E-13;
		break;
		
		case 449 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 450 :
		result = 1E-13;
		break;
		
		case 451 :
		result = 1E-13;
		break;
		
		case 452 :
		result = 1E-13;
		break;
		
		case 453 :
		result = 1E-13;
		break;
		
		case 454 :
		result = 1E-13;
		break;
		
		case 455 :
		result = 1E-13;
		break;
		
		case 456 :
		result = 1E-13;
		break;
		
		case 457 :
		result = 1E-13;
		break;
		
		case 458 :
		result = 1E-13;
		break;
		
		case 459 :
		result = 1E-13;
		break;
		
		case 460 :
		result = 1E-13;
		break;
		
		case 461 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 462 :
		result = 1E-13;
		break;
		
		case 463 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 464 :
		result = 1E-13;
		break;
		
		case 465 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 466 :
		result = 1E-13;
		break;
		
		case 467 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 468 :
		result = 1E-13;
		break;
		
		case 469 :
		result = 1E-13;
		break;
		
		case 470 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 471 :
		result = 1E-13;
		break;
		
		case 472 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 473 :
		result = 1E-13;
		break;
		
		case 474 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 475 :
		result = 1E-13;
		break;
		
		case 476 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 477 :
		result = 1E-13;
		break;
		
		case 478 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 479 :
		result = 1E-13;
		break;
		
		case 480 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 481 :
		result = 1E-13;
		break;
		
		case 482 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 483 :
		result = 1E-13;
		break;
		
		case 484 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 485 :
		result = 1E-13;
		break;
		
		case 486 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 487 :
		result = 1E-13;
		break;
		
		case 488 :
		result = 1E-13;
		break;
		
		case 489 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 490 :
		result = 1E-13;
		break;
		
		case 491 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 492 :
		result = 1E-13;
		break;
		
		case 493 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 494 :
		result = 1E-13;
		break;
		
		case 495 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 496 :
		result = 1E-13;
		break;
		
		case 497 :
		result = 2E-13 * pow(300/Tg,0.5);
		break;
		
		case 498 :
		result = 1E-13;
		break;
		
		case 499 :
		result = 1E-13;
		break;
		
		case 500 :
		result = 1E-13;
		break;
		
		case 501 :
		result = 1E-13;
		break;
		
		case 502 :
		result = 1E-13;
		break;
		
		case 503 :
		result = 1E-13;
		break;
		
		case 504 :
		result = 1E-13;
		break;
		
		case 505 :
		result = 1E-13;
		break;
		
		case 506 :
		result = 1E-13;
		break;
		
		case 507 :
		result = 1E-13;
		break;
		
		case 508 :
		result = 1E-13;
		break;
		
		case 509 :
		result = 8.3E-46 * exp(500/Tg);
		break;
		
		case 510 :
		result = 2.1E-17 * exp(100/Tg);
		break;
		
		case 511 :
		result = 5.8E-18 * exp(220/Tg);
		break;
		
		case 512 :
		result = 9.1E-19;
		break;
		
		case 513 :
		result = 6E-19;
		break;
		
		case 514 :
		result = 7E-19;
		break;
		
		case 515 :
		result = 6.3E-45 * exp(140/Tg);
		break;
		
		case 516 :
		result = 1.5E-17 * exp(-3600/Tg);
		break;
		
		case 517 :
		result = 5E-22;
		break;
		
		case 518 :
		result = 7.5E-17;
		break;
		
		case 519 :
		result = 1.7E-17 * exp(-1000/Tg);
		break;
		
		case 520 :
		result = 5E-18 * exp(-1620/Tg);
		break;
		
		case 521 :
		result = 1.5E-17 * exp(-570/Tg);
		break;
		
		case 522 :
		result = 6E-17;
		break;
		
		case 523 :
		result = 4.5E-17;
		break;
		
		case 524 :
		result = 7E-19;
		break;
		
		case 525 :
		result = 1.5E-18 * pow(Tg/300,0.5);
		break;
		
		case 526 :
		result = 6E-18 * pow(Tg/300,0.5);
		break;
		
		case 527 :
		result = 2.2E-20;
		break;
		
		case 528 :
		result = 4E-16;
		break;
		
		case 529 :
		result = 8E-17;
		break;
		
		case 530 :
		result = 8E-17;
		break;
		
		case 531 :
		result = 1.3E-17;
		break;
		
		case 532 :
		result = 5E-18 * exp(-210/Tg);
		break;
		
		case 533 :
		result = 1E-18;
		break;
		
		case 534 :
		result = 5E-17;
		break;
		
		case 535 :
		result = 7E-18;
		break;
		
		case 536 :
		result = 2.3E-17;
		break;
		
		case 537 :
		result = 5E-20;
		break;
		
		case 538 :
		result = 5E-17;
		break;
		
		case 539 :
		result = 1.25E5; ///???
		break;
		
		case 540 :
		result = 2.4E-16;
		break;
		
		case 541 :
		result = 3E-16;
		break;
		
		case 542 :
		result = 2.5E-17;
		break;
		
		case 543 :
		result = 3.2E-47 * exp(900/Tg);
		break;
		
		case 544 :
		result = 3.4E-46 * pow(300/Tg,1.2);
		break;
		
		case 545 :
		result = 8E-18 * exp(-2060/Tg);
		break;
		
		case 546 :
		result = 1E-43 * pow(300/Tg,1.6);
		break;
		
		case 547 :
		result = 6.5E-18 * exp(120/Tg);
		break;
		
		case 548 :
		result = 9E-44 * pow(300/Tg,2);
		break;
		
		case 549 :
		result = 1.7E-17;
		break;
		
		case 550 :
		result = 1.62E-44;
		break;
		
		case 551 :
		result = 2.2E-17 * exp(-350/Tg);
		break;
		
		case 552 :
		result = 3.3E-17 * exp(-2950/Tg);
		break;
		
		case 553 :
		result = 8.3E-17 * exp(-500/Tg);
		break;
		
		case 554 :
		result = 5.99E-17;
		break;
		
		case 555 :
		result = 2E-17 * exp(-3000/Tg);
		break;
		
		case 556 :
		result = 6.4E-18 * exp(67/Tg);
		break;
		
		case 557 :
		result = 8E-18;
		break;
		
		case 558 :
		result = 1E-17;
		break;
		
		case 559 :
		result = 1E-18;
		break;
		
		case 560 :
		result = 1.2E-16;
		break;
		
		case 561 :
		result = 1.2E-16;
		break;
		
		case 562 :
		result = 1.8E-17 * exp(107/Tg);
		break;
		
		case 563 :
		result = 9E-49;
		break;
		
		case 564 :
		result = 4.4E-17;
		break;
		
		case 565 :
		result = 7.2E-17;
		break;
		
		case 566 :
		result = 4E-17;
		break;
		
		case 567 :
		result = 1.4E-16;
		break;
		
		case 568 :
		result = 1.1E-16;
		break;
		
		case 569 :
		result = 2.2E-16;
		break;
		
		case 570 :
		result = 3.8E-24 * exp(-205/Tg);
		break;
		
		case 571 :
		result = 5.2E-17 * exp(-2840/Tg);
		break;
		
		case 572 :
		result = 8E-27;
		break;
		
		case 573 :
		result = 2.5E-17;
		break;
		
		case 574 :
		result = 1.5E-24;
		break;
		
		case 575 :
		result = 1.8E-18 * exp(-1370/Tg);
		break;
		
		case 576 :
		result = 1.4E-19 * exp(-2470/Tg);
		break;
		
		case 577 :
		result = 3.92E-16 * exp(-11400/Tg);
		break;
		
		case 578 :
		result = 2.8E-17 * pow(Tg/300,0.75);
		break;
		
		case 579 :
		result = 1.6E-18 * exp(-1000/Tg);
		break;
		
		case 580 :
		result = 1.4E-20 * exp(-600/Tg);
		break;
		
		case 581 :
		result = 3.09E-46 * pow(300/Tg,7.7);
		break;
		
		case 582 :
		result = 1.8E-17 * exp(110/Tg);
		break;
		
		case 583 :
		result = 7.4E-43 * pow(300/Tg,2.4);
		break;
		
		case 584 :
		result = 1E-44 * exp(300/Tg);
		break;
		
		case 585 :
		result = 3.4E-18 * exp(270/Tg); 
		break;
		
		case 586 :
		result = 3.3E-19 * exp(-1000/Tg);
		break;
		
		case 587 :
		result = 1.17E-45 * pow(300/Tg,3.8);
		break;
		
		case 588 :
		result = 2.8E-42 * pow( (300/Tg), 3.5);
		break;
		
		case 589 :
		result = 2.3E-19 * exp(-1600/Tg);
		break;
		
		case 590 :
		result = 1.47E-16;
		break;
		
		case 591 :
		result = 2.2E-42 * pow(300/Tg,2.9);
		break;
		
		case 592 :
		result = 5E-18 * exp(-3000/Tg);
		break;
		
		case 593 :
		result = 5.8E-16 * exp(-750/Tg);
		break;
		
		case 594 :
		result = 2E-17;
		break;
		
		case 595 :
		result = 4.8E-18;
		break;
		
		case 596 :
		result = 9.2E-19;
		break;
		
		case 597 :
		result = 1.03E-16 * exp(-2628/Tg);
		break;
		
		case 598 :
		result = 1.09E-13 * exp(-4952/Tg);
		break;
		
		case 599 :
		result = 1E-9 * exp(-11000/Tg) * pow( (300/Tg), 3.5);
		break;
		
		case 600 :
		result = 5.4E-44 * pow(Tg/300,-1.8);
		break;
		
		case 601 :
		result = 1.8E-42 / Tg;
		break;
		
		case 602 :
		result = 6.1E-38 / pow(Tg,2);
		break;
		
		case 603 :
		result = 1.69E-17 * exp(-1800/Tg);
		break;
		
		case 604 :
		result = 2.8E-18 * exp(-1900/Tg);
		break;
		
		case 605 :
		result = 5.6E-18;
		break;
		
		case 606 :
		result = 2.4E-18;
		break;
		
		case 607 :
		result = 4.2E-16 * exp(-950/Tg);
		break;
		
		case 608 :
		result = 3E-17 * exp(-500/Tg);
		break;
		
		case 609 :
		result = 2E-17 * exp(-3700/Tg);
		break;
		
		case 610 :
		result = 1.39E-20 * pow(Tg/298,3.29) * exp(-3160/Tg);
		break;
		
		case 611 :
		result = 3.2E-17 * exp(-2600/Tg);
		break;
		
		case 612 :
		result = 8.8E-18 * exp(-503/Tg);
		break;
		
		case 613 :
		result = 6.9E-43 * pow(Tg/300,-0.8);
		break;
		
		case 614 :
		result = 4.8E-17 * exp(250/Tg);
		break;
		
		case 615 :
		result = 2.9E-18 * exp(-160/Tg);
		break;
		
		case 616 :
		result = 8E-17 * exp(-500/Tg);
		break;
		
		case 617 :
		result = 1.8E-17 * exp(-390/Tg);
		break;
		
		case 618 :
		result = 1.5E-20 * exp(650/Tg);
		break;
		
		case 619 :
		result = 2.2E-19 * exp(600/Tg);
		break;
		
		case 620 :
		result = 5.25E-18 * exp(-1510/Tg);
		break;
		
		case 621 :
		result = 1.66E-21;
		break;
		
		case 622 :
		result = 1.4E-21 * exp(-1600/Tg);
		break;
		
		case 623 :
		result = 1E-26;
		break;
		
		case 624 :
		result = 1.6E-23;
		break;

		// Peng reactions R625-R673: BOLSIG+ cross-section rates
		case 625: case 626: case 627: case 628: case 629: case 630:
		case 631: case 633: case 634: case 636:
		case 638: case 639: case 640: case 641: case 642:
		case 643: case 644: case 645: case 646: case 647:
		case 648: case 649: case 650: case 651: case 652: case 653:
		case 655: case 656: case 657: case 658: case 659:
		case 660: case 661: case 663: case 665: case 666:
		case 667: case 668: case 669: case 670: case 671:
		case 672: case 673:
		{
			double TeEv = kb * Te / qe;
			result = interpolateBolsigRate(j, TeEv);
		}
		break;

		// Peng reactions with no cross-section data: rate = 0
		case 632: case 635: case 637: case 654:
		case 662: case 664:
		result = 0;
		break;

		default :
		result =0;
	}
	return result;
}



void reaction::setReactionRate(int j, double Tgas, double Telectron)
{
	reactionRate = reactionRateFunction(j, Tgas, Telectron);
}



void reaction::setReactantAndProductSpecies(int j)
{
	switch(j)
    {
		case 1:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 28; //N(2_D)
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 17; //e
		break;
		
		case 2:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //n2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 29; //N2(A_3_Sigma)
		productSpeciesList[2] = 17; //e 
		break;
		
		case 3:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 30; //N2(B_3_Pi)
		productSpeciesList[2] = 17; //e
		break;
		
		case 4:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 2; //N2+
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 17; //e
		break;
		
		case 5:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 28; //N(2_D)
		productSpeciesList[2] = 17; //e
		break;
		
		case 6:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 1;  //N+
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 17; //e
		break;
		
		case 7:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 17; //e
		break;
		
		case 8:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 31; //O(1_D)
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 17; //e
		break;
		
		case 9:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 35; //O2(a_1_Delta)
		productSpeciesList[2] = 17; //e
		break;
		
		case 10:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 17; //e
		break;
		
		case 11:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 12:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 31; //O(1_D)
		productSpeciesList[2] = 17; //e
		break;
		
		case 13:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 17; //e
		break;
		
		case 14:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 17; //e
		break;
		
		case 15:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 14; //OH+
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 17; //e
		productSpeciesList[4] = 17; //e
		break;
		
		case 16:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 11; //H+
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 17; //e
		productSpeciesList[4] = 17; //e
		break;
		
		case 17:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 44; //H2
		productSpeciesList[3] = 17; //e
		productSpeciesList[4] = 17; //e
		break;
		
		case 18:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 12; //H2+
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 17; //e
		productSpeciesList[4] = 17; //e
		break;
		
		case 19:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 17; //e
		break;
		
		case 20:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 44; //H2
		productSpeciesList[2] = 31; //O(1_D)
		productSpeciesList[3] = 17; //e
		break;
		
		case 21:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 17; //e
		break;
		
		case 22:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 12; //H2+
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 17; //e
		break;
		
		case 23:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 43; //N2O5
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 40; //NO3
		productSpeciesList[3] = 17; //e
		productSpeciesList[4] = 17; //e
		break;
		
		case 24:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 1;  //N+
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 0;  //M
		break;
		
		case 25:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 2;  //N2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 33; //N
		break;
		
		case 26:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 2;  //N2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 28; //N(2_D) 
		productSpeciesList[2] = 33; //N
		break;
		
		case 27:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 2;  //N2+
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2s
		productSpeciesList[2] = 0;  //M
		break;
		
		case 28:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 3;  //N3+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 33; //N
		break;
		
		case 29:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 4;  //N4+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 51; //N2
		break;
		
		case 30:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 5;  //O+
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 0;  //M
		break;
		
		case 31:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 6;  //O2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		break;
		
		case 32:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 6;  //O2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 31; //O(1_D)
		break;
		
		case 33:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 6;  //O2+
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 0;  //M
		break;
		
		case 34:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 7;  //O4+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		break;
		
		case 35:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 9;  //N2O+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		break;
		
		case 36:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 8;  //NO+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 34; //O
		break;
		
		case 37:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 8;  //NO+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 28; //N(2_D)
		productSpeciesList[2] = 34; //O
		break;
		
		case 38:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 8;  //NO+
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 0;  //M
		break;
		
		case 39:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 10; //NO2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 34; //O
		break;
		
		case 40:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 12; //H2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 32; //H
		break;
		
		case 41:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 13; //H3+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 44; //H2
		break;
		
		case 42:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 13; //H3+
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;
		
		case 43:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 15; //H2O+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		break;
		
		case 44:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 15; //H2O+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 44; //H2
		break;
		
		case 45:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 15; //H2O+
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;
		
		case 46:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 16; //H3O+
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;
		
		case 47:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 1;  //N+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 17; //e
		break;
		
		case 48:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 2;  //N2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 17; //e
		break;
		
		case 49:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 5;  //O+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 17; //e
		break;
		
		case 50:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e 
		reactantSpeciesList[3] = 6;  //O2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 17; //e
		break;
		
		case 51:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 8;  //NO+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 17; //e
		break;
		
		case 52:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 11; //H+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 17; //e
		break;
		
		case 53:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 12; //H2+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 44; //H2
		productSpeciesList[2] = 17; //e
		break;
		
		case 54:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 14; //OH+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 17; //e
		break;
		
		case 55:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 17; //e
		reactantSpeciesList[3] = 15; //H2O+
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 53; //H2O
		productSpeciesList[2] = 17; //e
		break;
		
		case 56:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 34; //O
		reactantSpeciesList[3] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 52; //O2
		break;
		
		case 57:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 34; //O
		reactantSpeciesList[3] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 34; //O
		break;
		
		case 58:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 52; //O2
		break;
		
		case 59:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 51; //N2
		break;
		
		case 60:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 34; //O
		break;
		
		case 61:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 1;
		productSpeciesList[1] = 19; //O2-
		break;
		
		case 62:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 52; //O2
		break;
		
		case 63:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 34; //O
		break;
		
		case 64:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 36; //O3
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 0;  //M
		break;
		
		case 65:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 51; //N2
		break;
		
		case 66:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 37; //NO
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 22; //NO-
		productSpeciesList[2] = 0;  //M
		break;
		
		case 67:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 37; //NO
		break;
		
		case 68:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 39; //NO2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 0;  //M
		break;
		
		case 69:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 40; //NO3
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 0;  //M
		break;
		
		case 70:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 26; //H-
		productSpeciesList[2] = 45; //OH
		break;
		
		case 71:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 44; //H2
		break;
		
		case 72:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 32; //H
		break;
		
		case 73:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 50; //HNO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 45; //OH
		break;
		
		case 74:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 17; //e
		break;
		
		case 75:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 17; //e
		break;
		
		case 76:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 29; //N2(A_3_Sigma)
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 17; //e
		break;
		
		case 77:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 30; //N2(B_3_Pi)
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 17; //e
		break;
		
		case 78:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 17; //e
		break;
		
		case 79:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 17; //e
		break;
		
		case 80:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 35; //O2(a_1_Delta)
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 17; //e
		break;
		
		case 81:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 82:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 17; //e
		break;
		
		case 83:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 53; //H2O
		productSpeciesList[2] = 17; //e
		break;
		
		case 84:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 17; //e
		break;
		
		case 85:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 86:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 29; //N2(A_3_Sigma)
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 87:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 30; //N2(B_3_Pi)
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 88:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 17; //e
		break;
		
		case 89:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 90:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 91:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 35; //O2(a_1_Delta)
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 92:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 17; //e
		break;
		
		case 93:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 94:
		reactantSpeciesList[0] = 2; 
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 95:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		productSpeciesList[4] = 17; //e
		break;
		
		case 96:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 0;  //M
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 0;  //M
		productSpeciesList[3] = 17; //e
		break;
		
		case 97:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 37; //NO
		productSpeciesList[3] = 17; //e
		break;
		
		case 98:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 38; //N2O
		productSpeciesList[3] = 17; //e
		break;
		
		case 99:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 44; //H2
		productSpeciesList[3] = 17; //e
		break;
		
		case 100:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 17; //e
		break;
		
		case 101:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 102:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 25; //NO3-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 17; //e
		break;
		
		case 103:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 25; //NO3-
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 36; //O3
		productSpeciesList[3] = 17; //e
		break;
		
		case 104:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 26; //H-
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 17; //e
		break;
		
		case 105:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 26; //H-
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 44; //H2
		productSpeciesList[2] = 17; //e
		break;
		
		case 106:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 27; //OH-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 17; //e
		break;
		
		case 107:    
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 27; //OH-
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 53; //H2O
		productSpeciesList[2] = 17; //e
	    break;
		
		case 108:  
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 33; //N
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 2;  //N2+
		productSpeciesList[2] = 0;  //M
		break;
		
		case 109: 
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 51; //N2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 3;  //N3+
		productSpeciesList[2] = 0;  //M
		break;
		
		case 110:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 51; //N2
		break;
		
		case 111:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 33; //N
	    break;
		
		case 112:  
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 2;  //N2+
		productSpeciesList[2] = 34; //O
		break;
		
		case 113:   
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 51; //N2
		break;
		
		case 114:    
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 33; //N
		break;
		
		case 115: 
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 37; //NO
		break;
		
		case 116:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 33; //N
		break;
		
		case 117:       
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 34; //O
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 0;  //M
		break;
		
		case 118: 
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 34; //O
		break;
		
		case 119:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+ 
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 37; //NO
		break;
		
		case 120:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 33; //N
		break;
		
		case 121:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 52; //O2
		break;
		
		case 122:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 14; //OH+
		productSpeciesList[2] = 33; //N
		break;
		
		case 123:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 32; //H
		break;
		
		case 124:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 33; //N
		break;
		
		case 125:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 44; //H2
		break;
		
		case 126:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 1;  //N+
		productSpeciesList[2] = 51; //N2
		break; 
		
		case 127:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 33; //N
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 3;  //N3+
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 128:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 51; //N2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 4;  //N4+
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 129:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 29; //N2(A_3_Sigma)
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 3;  //N3+
		productSpeciesList[2] = 33; //N
		break;   
		
		case 130:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 9;  //N2O+
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 131:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 132:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 133:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 134:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 135:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 33; //N
		break;   
		
		case 136:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 28; //N(2_D)
		break;   
		
		case 137:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 138:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 139:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 140:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 51; //N2
		break;
		
		case 141: 
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 2;  //N2+
		productSpeciesList[2] = 51; //N2
		break;
		
		case 142:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 51; //N2
		break;
		
		case 143:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 51; //N2
		break;
		
		case 144:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 51; //N2
		break;
		
		case 145:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;
		
		case 146:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 9;  //N2O+
		productSpeciesList[2] = 51; //N2
		break;
		
		case 147:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 148:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 37; //NO
		productSpeciesList[3] = 51; //N2
		break;
		
		case 149:   
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 51; //N2
		break;
		
		case 150:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 2;  //N2+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 151:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 9;  //N2O+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 152:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 1;  //N+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 153:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 154:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 155:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 38; //N2O
		productSpeciesList[3] = 51; //N2
		break;
		
		case 156:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 157:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 158:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 12; //H2+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 159:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;
		
		case 160:       
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 33; //N
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 0;  //M
		break;
		
		case 161: 
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 51; //N2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 0;  //M
		break;
		
		case 162:       
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 34; //O
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 0;  //M
		break;
		
		case 163:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 34; //O
		break;
		
		case 164:       
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 52; //O2
		break;
		
		case 165:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 28; //N(2_D)
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 1;  //N+
		productSpeciesList[2] = 34; //O
		break;
		
		case 166:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 34; //O
		break;
		
		case 167:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 33; //N
		break;   
		
		case 168:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 9;  //N2O+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 169:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5; //O+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 170:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 171:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 172:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 173:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 11; //H+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 174:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 14; //OH+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 175:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 14; //OH+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 176:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 177:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 178:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 7;  //O4+
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 179:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 180:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 181:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 43; //N2O5
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 40; //NO3
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 182:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 183:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 184:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 36; //O3
		break;   
		
		case 185:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 36; //O3
		break;   
		
		case 186:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 187:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 188:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 10; //NO2+ 
		productSpeciesList[2] = 52; //O2 
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 189:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 190:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 35; //O2(a_1_Delta)
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 34; //O
		break;   
		
		case 191:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 34; //O
		break;   
		
		case 192:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 22; //NO-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 193:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 23; //N2O-
		productSpeciesList[2] = 34; //O
		break;   
		
		case 194:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 37; //NO
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 195:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 34; //O
		break;   
		
		case 196:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 34; //O
		break;   
		
		case 197:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 32; //H
		break;   
		
		case 198:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 18; //O-
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 199:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 200:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 21; //O4-
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 201:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 202:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 203:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 204:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 205:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 19; //O2-
		reactantSpeciesList[2] = 50; //HNO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 46; //HO2
		break;   
		
		case 206:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 207:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 208:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 34; //O
		break;   
		
		case 209:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 210:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 36; //O3
		break;   
		
		case 211:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3- 
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 36; //O3
		break;   
		
		case 212:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 20; //O3-
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 213:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 21; //O4-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 214:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 21; //O4-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 215:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 21; //O4-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 216:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 21; //O4-
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 217:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 218:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 37; //NO
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 219:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 220:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 221:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 222:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 223:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 43; //N2O5
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 39; //NO2
		productSpeciesList[3] = 39; //NO2
		break;   
		
		case 224:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 225:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 33; //N
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 9;  //N2O+
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 226:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 227:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 228:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 229:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 19; //O2-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 230:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 231:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 22; //NO-
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 18; //O-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 232:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 233:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 43; //N2O5
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 40; //NO3
		productSpeciesList[3] = 37; //NO
		break;   
		
		case 234:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 22; //NO-
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 235:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 236:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 237:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 238:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 239:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 240:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 24; //NO2-
		reactantSpeciesList[2] = 50; //HNO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 25; //NO3-
		productSpeciesList[2] = 49; //HNO2
		break;   
		
		case 241:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 25; //NO3-
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 242:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 5;  //O+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 243:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 244:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 245:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 44; //H2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 13; //H3+
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 246:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+ 
		productSpeciesList[2] = 32; //H
		break;   
		
		case 247:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 26; //H-
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 248:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 26; //H-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 32; //H
		break;   
		
		case 249:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 26; //H-
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 27; //OH-
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 250:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 251:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 11; //H+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 252:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 13; //H3+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 253:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 16; //H3O+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 254:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 255:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 14; //OH+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 256:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 16; //H3O+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 257:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 258:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 259:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 260:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 46; //HO2
		break;   
		
		case 261:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 9;  //N2O+
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 262:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 263:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 264:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 15; //H2O+
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 265:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 16; //H3O+
		productSpeciesList[2] = 34; //O
		break;   
		
		case 266:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 27; //OH-
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 20; //O3-
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 267:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 27; //OH-
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 24; //NO2-
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 268:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 269:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 270:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 6;  //O2+
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 271:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 272:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 273:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 16; //H3O+
		productSpeciesList[2] = 32; //H
		break;   
		
		case 274:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 16; //H3O+
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 275:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 8;  //NO+
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 53; //H2O
		break;   
		
		case 276:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 43; //N2O5
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 10; //NO2+
		productSpeciesList[2] = 50; //HNO3
		productSpeciesList[3] = 53; //H20
		break;   
		
		case 277:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		break;   
		
		case 278:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 34; //O
		break;   
		
		case 279:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 34; //O
		break;   
		
		case 280:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 281:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 34; //O
		break;   
		
		case 282:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 34; //O
		break;   
		
		case 283:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 34; //O
		break;   
		
		case 284:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 34; //O
		break;   
		
		case 285:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 32; //H
		break;   
		
		case 286:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 5;  //O+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 287:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 288:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 289:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 290:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 291:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 292:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 293:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 294:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 295:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 296:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 297:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 298:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 299:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 300:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 301:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 302:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 303:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 304:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 305:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 6;  //O2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 306:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 307:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 308:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 309:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		productSpeciesList[4] = 52; //O2
		break;   
		
		case 310:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 311:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 312:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 313:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 314:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 315:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 7;  //O4+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 316:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 33; //N
		break;   
		
		case 317:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 33; //N
		break;   
		
		case 318:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 33; //N
		break;   
		
		case 319:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 33; //N
		break;   
		
		case 320:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 33; //N
		break;   
		
		case 321:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 33; //N
		break;   
		
		case 322:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 33; //N
		break;   
		
		case 323:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 324:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 32; //H
		break;   
		
		case 325:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 1;  //N+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 326:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 327:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 328:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 329:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 330:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 331:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 332:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 333:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 334:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 335:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 336:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 337:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 338:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 339:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 340:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 341:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 342:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 343:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 344:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 2;  //N2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;   
		
		case 345:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 346:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 347:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 348:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 349:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 350:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 351:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 352:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 353:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 354:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 355:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 356:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 357:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 358:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 359:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 38; //N2O
		break;   
		
		case 360:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 361:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 362:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 38; //N2O
		break;   
		
		case 363:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 9;  //N2O+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 34; //O
		break;   
		
		case 364:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 365:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 366:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 367:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 368:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 369:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 370:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 371:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 33; //N
		productSpeciesList[4] = 51; //N2
		break;   
		
		case 372:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 373:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 3;  //N3+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 33; //N
		break;   
		
		case 374:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 375:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 376:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 377:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 378:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 379:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 380:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 381:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 51; //N2
		productSpeciesList[4] = 51; //N2
		break;   
		
		case 382:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 383:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 4;  //N4+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 384:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 385:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 386:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 387:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 388:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 389:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 390:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+ 
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 391:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 392:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 393:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 394:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 395:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 396:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 397:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 398:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 37; //NO
		break;   
		
		case 399:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 400:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 401:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 402:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 8;  //NO+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 34; //O
		break;   
		
		case 403:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 404:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 405:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 406:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 407:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 408:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 409:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 410:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 411:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 412:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 413:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 414:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 415:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 416:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 417:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 39; //NO2
		break;   
		
		case 418:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 419:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 420:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 421:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 10; //NO2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 422:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 32; //H
		break;   
		
		case 423:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 32; //H
		break;   
		
		case 424:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 32; //H
		break;   
		
		case 425:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 32; //H
		break;   
		
		case 426:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 32; //H
		break;   
		
		case 427:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 32; //H
		break;   
		
		case 428:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 32; //H
		break;   
		
		case 429:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 32; //H
		break;   
		
		case 430:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 32; //H
		break;   
		
		case 431:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 11; //H+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		break;   
		
		case 432:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 433:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 434:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 435:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 436:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 437:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 438:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 439:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 440:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 441:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 442:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 443:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 444:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 445:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 446:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 447:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 448:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 449:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 450:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 12; //H2+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 32; //H
		break;   
		
		case 451:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 452:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 453:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 454:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 455:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 456:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 457:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 458:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 32; //H
		productSpeciesList[4] = 44; //H2
		break;   
		
		case 459:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 460:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 13; //H3+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 44; //H2
		break;   
		
		case 461:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 462:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 463:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 464:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 465:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 466:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 467:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 468:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 469:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 470:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 471:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 472:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 473:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 474:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 475:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 45; //OH
		break;   
		
		case 476:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 477:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 478:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 479:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 14; //OH+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 480:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 481:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 482:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 483:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 484:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 485:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 486:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 487:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 488:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 489:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 490:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 491:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 492:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 493:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 494:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 53; //H2O
		break;   
		
		case 495:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 496:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 26; //H- 
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 497:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 53; //H2O
		break;   
		
		case 498:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 15; //H2O+
		reactantSpeciesList[2] = 27; //OH- 
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 32; //H
		break;   
		
		case 499:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 23; //N2O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 500:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 22; //NO-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 501:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 24; //NO2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 502:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 25; //NO3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 503:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 18; //O-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 504:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 19; //O2-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 505:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 20; //O3-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 506:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 21; //O4-
		productSpeciesList[0] = 4;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 53; //H2O
		productSpeciesList[4] = 32; //H
		break;   
		
		case 507:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 26; //H-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 508:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 16; //H3O+
		reactantSpeciesList[2] = 27; //OH-
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 53; //H2O
		productSpeciesList[3] = 32; //H
		break;   
		
		case 509:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 33; //N
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 510:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		break;   
		
		case 511:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 34; //O
		break;   
		
		case 512:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 513:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 514:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 515:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 34; //O
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 516:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 34; //O
		break;   
		
		case 517:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 518:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 519:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 33; //N
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 520:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 521:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 522:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 1;
		productSpeciesList[1] = 38; //N2O
		break;   
		
		case 523:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 524:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 34; //O
		break;   
		
		case 525:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 34; //O
		break;   
		
		case 526:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 28; //N(2_D)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 31; //O(1_D)
		break;   
		
		case 527:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 528:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 29; //N2(A_3_Sigma)
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 30; //N2(B_3_Pi)
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 529:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 530:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 531:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 37; //NO
		productSpeciesList[3] = 34; //O
		break;   
		
		case 532:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 533:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 35; //O2(a_1_Delta)
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 534:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 33; //N
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 33; //N
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 535:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 28; //N(2_D)
		break;   
		
		case 536:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 31; //O(1_D)
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 537:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 29; //N2(A_3_Sigma)
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 51; //N2
		break;   
		
		case 538:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 30; //N2(B_3_Pi)
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 29; //N2(A_3_Sigma)
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 539:
		reactantSpeciesList[0] = 1;
		reactantSpeciesList[1] = 30; //N2(B_3_Pi)
		productSpeciesList[0] = 1;
		productSpeciesList[1] = 29; //N2(A_3_Sigma)
		break;   
		
		case 540:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 30; //N2(B_3_Pi)
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 29; //N2(A_3_Sigma)
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 541:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 30; //N2(B_3_Pi)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;   
		
		case 542:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 30; //N2(B_3_Pi)
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 29; //N2(A_3_Sigma)
		productSpeciesList[2] = 44; //H2
		break;   
		
		case 543:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 34; //O
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 544:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 36; //O3
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 545:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 546:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 37; //NO
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 547:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O 
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 548:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 39; //NO2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 549:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 550:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 32; //H
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 551:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 552:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 47; //H2O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 46; //HO2
		break;   
		
		case 553:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 554:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 48; //HNO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 555:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 34; //O
		reactantSpeciesList[2] = 49; //HNO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 556:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 557:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 34; //O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		break;   
		
		case 558:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 35; //O2(a_1_Delta)
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 559:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 35; //O2(a_1_Delta)
		break;   
		
		case 560:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 36;
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 52; //O2
		break;  
		
		case 561:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 562:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 563:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 51; //N2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 564:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 51; //N2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 565:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 38; //N2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 566:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 567:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 568:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 44; //H2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 32; //H
		break;   
		
		case 569:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 31; //O(1_D)
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 45; //OH
		break;   
		
		case 570:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 35; //O2(a_1_Delta)
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 571:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 35; //O2(a_1_Delta)
		reactantSpeciesList[2] = 36; //O3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 572:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 35; //O2(a_1_Delta)
		reactantSpeciesList[2] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 573:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 35; //O2(a_1_Delta)
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 37; //NOs
		break;   
		
		case 574:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 35; //O2(a_1_Delta)
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 51; //N2
		break;   
		
		case 575:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 36; //O3
		reactantSpeciesList[2] = 37; //NO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 576:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 36; //O3
		reactantSpeciesList[2] = 39; //NO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 577:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 36; //O3
		reactantSpeciesList[2] = 0;  //M
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 0;  //M
		break;   
		
		case 578:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 36; //O3
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 579:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 36; //O3
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 580:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 36; //O3
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 52; //O2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 581:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 37; //NO
		reactantSpeciesList[2] = 39; //NO2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 41; //N2O3
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 582:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 37; //NO
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 583:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 37; //NO
		reactantSpeciesList[2] = 45; //OH
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 49; //HNO2
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 584:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 37; //NO
		reactantSpeciesList[2] = 32; //H
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 48; //HNO
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 585:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 37; //NO
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 586:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 37; //NO
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 48; //HNO
		break;   
		
		case 587:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 39; //NO2
		reactantSpeciesList[2] = 39; //NO2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 42; //N2O4
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 588:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 39; //NO2
		reactantSpeciesList[2] = 40; //NO3
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 43; //N2O5
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 589:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 39; //NO2
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 37; //NO
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 590:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 39; //NO2
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 37; //NO
		break;   
		
		case 591:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 39; //NO2
		reactantSpeciesList[2] = 45; //OH
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 50; //HNO3
		productSpeciesList[2] = 0;  //M
		break;   
		
		case 592:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 40; //NO3
		reactantSpeciesList[2] = 40; //NO3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 39; //NO2
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 593:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 40; //NO3
		reactantSpeciesList[2] = 32; //H
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 594:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 40; //NO3
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 39; //NO2
		break;   
		
		case 595:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 40; //NO3
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 45; //OH
		productSpeciesList[3] = 52; //O2
		break;   
		
		case 596:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 40; //NO3
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 50; //HNO3
		productSpeciesList[2] = 52; //O2
		break;   
		
		case 597:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 41; //N2O3
		reactantSpeciesList[2] = 0;  //M
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO  
		productSpeciesList[2] = 39; //NO2
		productSpeciesList[3] = 0;  //M
		break;   
		
		case 598:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 42; //N2O4
		reactantSpeciesList[2] = 0;  //M
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 39; //NO2
		productSpeciesList[3] = 0;  //M
		break;   
		
		case 599:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 43; //N2O5
		reactantSpeciesList[2] = 0;  //M
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 40; //NO3
		productSpeciesList[3] = 0;  //M
		break; 
		
		case 600:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 52; //O2
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 0;  //M
		break;
		
		case 601:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 32; //H
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 44; //H2
		productSpeciesList[2] = 0;  //M
		break;
		
		case 602:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 45; //OH
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 53; //H2O
		productSpeciesList[2] = 0;  //M
		break;
		
		case 603:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 47; //H2O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 604:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 47; //H2O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 44; //H2
		break;
		
		case 605:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 44; //H2
		productSpeciesList[2] = 52; //O2
		break;
		
		case 606:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 607:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 45; //OH
		break;
		
		case 608:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 48; //HNO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 44; //H2
		break;
		
		case 609:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 49; //HNO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 44; //H2
		break;
		
		case 610:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 32; //H
		reactantSpeciesList[2] = 50; //HNO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 611:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 44; //H2
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 32; //H
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 612:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 45; //OH
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 613:
		reactantSpeciesList[0] = 3;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 45; //OH
		reactantSpeciesList[3] = 0;  //M
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 47; //H2O2
		productSpeciesList[2] = 0;  //M
		break;
		
		case 614:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 52; //O2
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 615:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 47; //H2O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 46; //HO2
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 616:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 48; //HNO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 617:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 49; //HNO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 618:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 45; //OH
		reactantSpeciesList[2] = 50; //HNO3
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 40; //NO3
		productSpeciesList[2] = 53; //H2O
		break;
		
		case 619:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 46; //HO2
		reactantSpeciesList[2] = 46; //HO2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 47; //H2O2
		productSpeciesList[2] = 52; //O2
		break;
		
		case 620:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 48; //HNO
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 46; //HO2
	    break;
		
		case 621:  
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 48; //HNO
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 45; //OH
		break;
		
		case 622:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 48; //HNO
		reactantSpeciesList[2] = 48; //HNO
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 38; //N2O
		productSpeciesList[2] = 53; //H2O
	    break;
		
		case 623:  
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 49; //HNO2
		reactantSpeciesList[2] = 49; //HNO2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 37; //NO
		productSpeciesList[2] = 39; //NO2
		productSpeciesList[3] = 53; //H2O
		break;
		
		case 624:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 49; //HNO2
		reactantSpeciesList[2] = 50; //HNO3
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 39; //NO2
		productSpeciesList[2] = 39; //NO2
		productSpeciesList[3] = 53; //H2O
		break;

		// === Peng reactions R625-R673 ===

		// R625: e + H2O -> H2 + O-  (attachment)
		case 625:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 44; //H2
		productSpeciesList[2] = 18; //O-
		break;

		// R626: e + H2O -> OH + H-  (attachment)
		case 626:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 45; //OH
		productSpeciesList[2] = 26; //H-
		break;

		// R627-R630: e + H2O -> e + H2O  (elastic/vibrational, energy-loss-only)
		case 627: case 628: case 629: case 630:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 53; //H2O
		break;

		// R631: e + H2O -> e + H + OH  (dissociation)
		case 631:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 32; //H
		productSpeciesList[3] = 45; //OH
		break;

		// R632: e + H2O -> e + H2O  (inverse, rate=0, energy-loss-only)
		case 632:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 53; //H2O
		break;

		// R633: e + H2O -> e + e + H2O+  (ionization)
		case 633:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 53; //H2O
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 15; //H2O+
		break;

		// R634-R644: e + N2 -> e + N2  (elastic/vibrational, energy-loss-only)
		case 634: case 635: case 636: case 637:
		case 638: case 639: case 640: case 641: case 642:
		case 643: case 644:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 51; //N2
		break;

		// R645: e + N2 -> e + N2(A)  (electronic exc. A3Sigma triplet)
		case 645:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 29; //N2(A)
		break;

		// R646: e + N2 -> e + N2(A)  (electronic exc. A3Sigma v1-4)
		case 646:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 29; //N2(A)
		break;

		// R647: e + N2 -> e + N2(B)  (electronic exc. B3Pi)
		case 647:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 30; //N2(B)
		break;

		// R648: e + N2 -> e + N2  (electronic exc. W3Delta, untracked)
		case 648:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 51; //N2
		break;

		// R649: e + N2 -> e + N2(A)  (electronic exc. A'1Sigma mapped to A)
		case 649:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 29; //N2(A)
		break;

		// R650-R656: e + N2 -> e + N2  (electronic exc. to untracked states)
		case 650: case 651: case 652: case 653:
		case 654: case 655: case 656:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 51; //N2
		break;

		// R657: e + N2 -> e + N + N  (dissociation)
		case 657:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 33; //N
		productSpeciesList[3] = 33; //N
		break;

		// R658: e + N2 -> e + e + N2+  (ionization)
		case 658:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 2;  //N2+
		break;

		// R659: e + N2 -> e + e + N2+  (ionization, second channel)
		case 659:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 51; //N2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 2;  //N2+
		break;

		// R660: e + O2 -> O + O-  (attachment)
		case 660:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 34; //O
		productSpeciesList[2] = 18; //O-
		break;

		// R661-R665: e + O2 -> e + O2  (elastic/vibrational, energy-loss-only)
		case 661: case 662: case 663: case 664: case 665:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 52; //O2
		break;

		// R666: e + O2 -> e + O2(a1D)  (electronic exc.)
		case 666:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 35; //O2(a1D)
		break;

		// R667-R668: e + O2 -> e + O2  (electronic exc. to untracked states)
		case 667: case 668:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 52; //O2
		break;

		// R669: e + O2 -> e + O + O  (dissociation)
		case 669:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 34; //O
		break;

		// R670: e + O2 -> e + O + O(1D)  (dissociation)
		case 670:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 34; //O
		productSpeciesList[3] = 31; //O(1D)
		break;

		// R671-R672: e + O2 -> e + O2  (electronic exc. to untracked states)
		case 671: case 672:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 2;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 52; //O2
		break;

		// R673: e + O2 -> e + e + O2+  (ionization)
		case 673:
		reactantSpeciesList[0] = 2;
		reactantSpeciesList[1] = 17; //e
		reactantSpeciesList[2] = 52; //O2
		productSpeciesList[0] = 3;
		productSpeciesList[1] = 17; //e
		productSpeciesList[2] = 17; //e
		productSpeciesList[3] = 6;  //O2+
		break;

		default:
		;
	}
}



int reaction::returnNumberOfReactants(void)
{ return reactantSpeciesList[0]; }

int reaction::returnNumberOfProducts(void)  
{ return productSpeciesList[0]; }

int reaction::returnReactant(int itemIndex)
{ return reactantSpeciesList[itemIndex]; }

int reaction::returnProduct(int itemIndex)
{ return productSpeciesList[itemIndex]; }

double reaction::returnReactionRate(void)
{ return reactionRate; }


