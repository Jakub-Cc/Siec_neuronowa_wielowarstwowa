#include "Siec.h"
#include <infstr.h>
#include <fstream>
#include <iostream>





int main()
{
	Tools::static_init();

	Tools::static_init();

	int ilosc_test = 0;
	int dlugosc_wekt = 0;
	double ** test = new double *[1];
	double ** test_etykiety = new double*[1];

	int ilosc_ucz = 0;
	double ** ucz = new double *[1];
	double ** ucz_etykiety = new double *[1];


	clock_t begin = clock();
	Tools::read_mnist_data("MNIST/t10k-images.idx3-ubyte", ilosc_test, dlugosc_wekt, test);
	Tools::read_mnist_labels("MNIST/t10k-labels.idx1-ubyte", ilosc_test, test_etykiety);
	Tools::read_mnist_data("MNIST/train-images.idx3-ubyte", ilosc_ucz, dlugosc_wekt, ucz);
	Tools::read_mnist_labels("MNIST/train-labels.idx1-ubyte", ilosc_ucz, ucz_etykiety);
	clock_t end = clock();

	std::cout << "Wczytano ciagi MNIST - Czas ladowania : " << double(end - begin) / CLOCKS_PER_SEC << "\n";

	//unsigned int fp_control_state = _controlfp ( _EM_INEXACT, _MCW_EM );


	Siec siec(dlugosc_wekt, 10, 100, 0.25);
	siec.ustaw_wagi(-0.5, 0.5);
	//Tools::zapis_do_pliku ( "siec0.txt", siec );

	Tools::odczyt_z_pliku ( "Wyniki/siec3reg.txt", siec );
	siec.zmien_ilosc_neuronow_war_wyj ( 10 );
	siec.ustaw_wagi_war_wyj(-0.5, 0.5);

	siec.jaka_f_akt_ukr = 0;
	siec.jaka_f_akt_wyj = 1;

	siec.ciagi_uczace = ucz;
	siec.ciagi_uczace_war_oczekiwana = ucz_etykiety;
	siec.ilosc_ciagow_uczacych = ilosc_ucz;
	siec.ciagi_walidacyjne = test;
	siec.ciagi_walidacyjne_war_oczekiwana = test_etykiety;
	siec.ilosc_ciagow_walidacyjnych = ilosc_test;

	std::ofstream MyExcelFile2;
	siec.reg_L2 = 1;
	siec.prawd_wylaczenia = 0;

	begin = clock();

	siec.epoki_uczenie ( 50, 0 );

	//siec.epoka ();
	end = clock();

	std::cout << "Czas : " << double(end - begin) / CLOCKS_PER_SEC << "\n";
	Tools::zapis_do_pliku("Wyniki/siec8.txt", siec);
	siec.save_to_csv("Wyniki/wagi.csv");
	Tools::wagi_to_png("Wyniki/test.png", siec.l, 28, 28, siec.wagi_war_ukr);
	//Tools::wagi_to_png_transp("Wyniki/test2.png", siec.l, 28, 28, siec.wagi_war_wyj);

	std::cout << "TEST: " << siec.test(test, test_etykiety, ilosc_test) << '\n';
	std::cout << "TEST2: " << siec.test2(test, test_etykiety, ilosc_test) << '\n';

	//MyExcelFile2.open ( "test_epok23.csv" );
	//MyExcelFile2 << "epoka;suma bledow;blad srednio kwad walidacji;\n";

	//siec.epoki_uczenie (2, 0 );
	//for ( int ab = 0; ab < 8; ab++ )
	{
		//	MyExcelFile2 << ab << ";" << siec.epoka () << ";" << siec.test2 ( test, test, ilosc_test ) << ";\n";
		//	Tools::zapis_do_pliku ( "S/siec" + std::to_string ( ab ) + ".txt", siec );
	}
	//MyExcelFile2.close ();
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";
	//std::cout << siec.epoka () << " " << siec.test2 ( test, test, ilosc_test ) << "\n";


	//Tools::zapis_do_pliku ( "siec.txt", siec );


	//zrzut wag do pliku
	//siec.dawaj_skal_wagi_ukr("wagi2.csv");
	//siec.dawaj_skal_wag_wyj("wagi3.csv");

	//maly test dzialania
	if (1)
	{
		//std::cout << siec.test2(test, test, ilosc_test) << "\n";
		std::ofstream MyExcelFile;
		MyExcelFile.open("Wyniki/test.csv");
		for (int k = 10; k < 15; k++)
		{

			for (int i = 0; i < 28; i++)
			{
				for (int j = 0; j < 28; j++)
				{
					MyExcelFile << test[k][i * 28 + j] << ";";
				}
				MyExcelFile << "\n";
			}
			MyExcelFile << "\n";

			std::streambuf* cout_sbuf = std::cout.rdbuf();
			std::cout.rdbuf(NULL);

			siec.policz_wyjscie(test[k], 2);

			std::cout.rdbuf(cout_sbuf);

			for (int j = 0; j < 10; j++)
			{
				MyExcelFile << test_etykiety[k][j] << ";";
			}

			MyExcelFile << "\n";

			for (int j = 0; j < 10; j++)
			{
				MyExcelFile << siec.wyj_war_wyj[j] << ";";

			}
			MyExcelFile << "\n";
			MyExcelFile << "\n";
		}
		MyExcelFile.close();
	}

	//big test
	if (1)
	{
		std::ofstream MyExcelFile;
		MyExcelFile.open("Wyniki/Big_test.csv");
		std::streambuf* cout_sbuf = std::cout.rdbuf();

		MyExcelFile << ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;=ŒREDNIA(AH3:AH" << ilosc_test + 2 << ");;=SUMA(AJ3:AJ" << ilosc_test + 2 << ")/" << ilosc_test << "\n";
		MyExcelFile << "wzorzec oczekiwany;;;;;;;;;;;wzorzec wyjsciowy;;;;;;;;;;;bledy srednio kwad na wyjsciach;;;;;;;;;;;suma bledow;;czy popr odp\n";

		for (int k = 0; k < ilosc_test; k++)
		{
			std::cout.rdbuf(NULL);
			siec.policz_wyjscie(test[k], 2);
			std::cout.rdbuf(cout_sbuf);

			for (int j = 0; j < 10; j++)
			{
				MyExcelFile << test_etykiety[k][j] << ";";
			}
			MyExcelFile << ";";

			for (int j = 0; j < 10; j++)
			{
				MyExcelFile << siec.wyj_war_wyj[j] << ";";
			}
			MyExcelFile << ";";

			for (int j = 0; j < 10; j++)
			{
				MyExcelFile << "=(" << (char)('a' + j) << k + 3 << "-" << (char)('a' + 11 + j) << k + 3 << ")^2;";
			}
			MyExcelFile << ";";
			MyExcelFile << "=0.5*SUMA(W" << k + 3 << ":AF" << k + 3 << ")";
			MyExcelFile << ";";
			MyExcelFile << ";";

			MyExcelFile << "= JE¯ELI(PODAJ.POZYCJÊ(MAX(A" << 3 + k << ":J" << 3 + k << "),A" << 3 + k << ":J" << 3 + k
				<< ",0)=PODAJ.POZYCJÊ(MAX(L" << 3 + k << ":U" << 3 + k << "),L" << 3 + k << ":U" << 3 + k << ", 0 ),1,0);";
			MyExcelFile << "\n";

		}

		MyExcelFile.close();
	}

	//heatmap pixeli
	if (0)
	{
		int ilosc[10];
		double *tab[10];
		double *mix=new double[dlugosc_wekt]();

		for (int i = 0; i < 10; i++)
		{
			tab[i] = new double[dlugosc_wekt]();
			ilosc[i] = 0;
		}
	
		for (int i = 0; i < ilosc_test; i++)
		{
			int cyfra = std::max_element(test_etykiety[i], test_etykiety[i] + 10)- test_etykiety[i];
			ilosc[cyfra]++;
			
			for (int j = 0; j < dlugosc_wekt; j++)
			{
				tab[cyfra][j] += test[i][j];
				mix[j] += test[i][j];
			}
		}

		std::ofstream MyExcelFile;
		MyExcelFile.open("Wyniki/heatmap.csv.csv");
		for (int i = 0; i < 10; i++)
		{
			MyExcelFile << ilosc[i] << "\n\n";
			for (int j = 0; j < 28; j++)
			{
				for (int k = 0; k < 28; k++)
				{
					MyExcelFile << tab[i][j * 28 + k] << ";";
				}
				MyExcelFile << "\n";
			}
			MyExcelFile << "\n";
		}
		for (int j = 0; j < 28; j++)
		{
			for (int k = 0; k < 28; k++)
			{
				MyExcelFile << mix[j * 28 + k] << ";";
			}
			MyExcelFile << "\n";
		}
		MyExcelFile << "\n";

		MyExcelFile.close();
	}


	system("pause");
	return 1;
}