#include "catch.hpp"
#include "config.h"
#include "vtk_writer.h"

using namespace learnSPH;

TEST_CASE("SimpleDam_Densities") {
	SECTION("SimpleDam_Densities") {
		std::stringstream outputFile;
		outputFile << SOURCE_DIR << "/res/densities/simpleDam";
		std::stringstream inputFileRoot;
		inputFileRoot << SOURCE_DIR << "/res/simple_dam/simpleDam";
		SECTION("SPH") {
			outputFile << "_SPH";
			inputFileRoot << "_SPH";
			SECTION("I") {
				outputFile << "_I_densities.csv";
				inputFileRoot << "_I";
			}
			SECTION("II") {
				outputFile << "_II_densities.csv";
				inputFileRoot << "_II";
			}
			SECTION("III") {
				outputFile << "_III_densities.csv";
				inputFileRoot << "_III";
			}
			SECTION("IV") {
				outputFile << "_IV_densities.csv";
				inputFileRoot << "_IV";
			}
			SECTION("V") {
				outputFile << "_V_densities.csv";
				inputFileRoot << "_V";
			}
			SECTION("VI") {
				outputFile << "_VI_densities.csv";
				inputFileRoot << "_VI";
			}
			SECTION("VII") {
				outputFile << "_VII_densities.csv";
				inputFileRoot << "_VII";
			}
			SECTION("VIII") {
				outputFile << "_VIII_densities.csv";
				inputFileRoot << "_VIII";
			}
		}

		SECTION("PBF") {
			outputFile << "_PBF";
			inputFileRoot << "_PBF";
			SECTION("I") {
				outputFile << "_I_densities.csv";
				inputFileRoot << "_I";
			}
			SECTION("II") {
				outputFile << "_II_densities.csv";
				inputFileRoot << "_II";
			}
			SECTION("III") {
				outputFile << "_III_densities.csv";
				inputFileRoot << "_III";
			}
			SECTION("IV") {
				outputFile << "_IV_densities.csv";
				inputFileRoot << "_IV";
			}
			SECTION("V") {
				outputFile << "_V_densities.csv";
				inputFileRoot << "_V";
			}
			SECTION("VI") {
				outputFile << "_VI_densities.csv";
				inputFileRoot << "_VI";
			}
			SECTION("VII") {
				outputFile << "_VII_densities.csv";
				inputFileRoot << "_VII";
			}
			SECTION("VIII") {
				outputFile << "_VIII_densities.csv";
				inputFileRoot << "_VIII";
			}
		}

		std::stringstream densitiesCSV;
		for (int i = 0; i <= 160; i++) {
			std::stringstream inputFile;
			inputFile << inputFileRoot.str() << std::to_string(i) << ".vtk";
			std::vector<double> densities;
			readDensitiesFromVTK(inputFile.str(), densities);
			double maxDensity = *std::max_element(densities.begin(), densities.end());
			densitiesCSV << std::to_string(maxDensity) << "\n";
		}

		std::ofstream outfile;
		outfile.open(outputFile.str());
		outfile << densitiesCSV.str();
		outfile.close();
	}
}