#include "ExpectedValuesForInternalCoordinatesTest.h"

//#include "../../Scon/scon_mathmatrix.h"

namespace ExpectedValuesForInternalCoordinates {

	std::unique_ptr<scon::mathmatrix<double>> ReadMatrixFiles::readNLinesOfFileWithMNumbers(std::size_t const n, std::size_t const m) {
		auto  nextMatrix = std::make_unique<scon::mathmatrix<double>>(n, m);

		for (auto i = 0u; i < n; ++i) {
			std::string line;
			std::getline(inputStream, line);
			std::stringstream lineStream(line);

			for (auto j = 0u; j < m; ++j) lineStream >> nextMatrix->operator()(i, j);
		}

		return nextMatrix;
	}

	/*std::unique_ptr<scon::mathmatrix<double>> expectedValuesForF() {
		return std::make_unique<scon::mathmatrix<double>>(std::initializer_list<std::initializer_list<double>>{ 
			std::initializer_list<double>{ 6.84606631495125, -4.57230316634422, -13.0165563000258, 3.74086935468525 },
			std::initializer_list<double>{ -4.57230316634422, 0.41134631529563, 9.81533983441827, -19.2451004432376 },
			std::initializer_list<double>{ -13.0165563000258, 9.81533983441827, -8.15722162019127, -10.0102206984924 },
			std::initializer_list<double>{ 3.74086935468525, -19.2451004432376, -10.0102206984924, 0.899808989944389 } });
	}

	std::unique_ptr<scon::mathmatrix<double>> expectedEigenValuesForF() {
		return std::make_unique<scon::mathmatrix<double>>(std::initializer_list<std::initializer_list<double>>{
			std::initializer_list<double>{ -18.7538392250112, },
			std::initializer_list<double> { -18.0459463050214, },
			std::initializer_list<double>{ 6.10328376463374, }, 
			std::initializer_list<double>{ 30.6965017653989, } });
	}

	std::unique_ptr<scon::mathmatrix<double>> expectedEigenVectorsForF() {
		return std::make_unique<scon::mathmatrix<double>>(std::initializer_list<std::initializer_list<double>>{ 
			std::initializer_list<double>{ -0.199380447172903, 0.342091987934248, -0.811125403595648, 0.430460321885915 },
			std::initializer_list<double>{ -0.549673535066269, -0.495167463863083, -0.371435099023426, -0.560984986245274 },
			std::initializer_list<double>{ -0.401611672479942, 0.783069660011025, 0.200535739330358, -0.430459509535228 },
			std::initializer_list<double>{ -0.70485069813454, -0.156793373880449, 0.404841900136933, 0.56098517550819 } });
	}

	std::unique_ptr<scon::mathmatrix<double>> expectedRotationDerivatives() {
		return std::make_unique<scon::mathmatrix<double>>(std::initializer_list<std::initializer_list<double>>{
			std::initializer_list<double>{ 0.0440694562731686, -0.0462415234932995, 0.0342224408171854 },
			std::initializer_list<double>{ 0.0250224091795194, 0.0582434877694512, 0.0342425808649325 },
			std::initializer_list<double>{ 0.023893908788188, 0.0013777688386713, 0.0231910960119916 },
			std::initializer_list<double>{ -0.00880081128293551, 0.09114820691289, -0.251807550239404 },
			std::initializer_list<double>{ -0.153659042727933, -0.15948765638584, -0.0393920943151345 },
			std::initializer_list<double>{ 0.0344559570944513, 0.0437853128846211, -0.0111303197288723 },
			std::initializer_list<double>{ -0.156235779332275, 0.0260477069160252, 0.291047355995862 },
			std::initializer_list<double>{ -0.0888487906434072, 0.154367319708192, 0.189852131893457 },
			std::initializer_list<double>{ -0.0846724302441519, -0.195572120300373, 0.121164438538402 },
			std::initializer_list<double>{ 0.0615615310599121, -0.0977936053539417, 0.14708852204346 },
			std::initializer_list<double>{ 0.344007690491845, 0.0300320610974034, -0.193800748922808 },
			std::initializer_list<double>{ -0.0481724531386894, 0.0998496802546184, -0.15619456025363 },
			std::initializer_list<double>{ 0.280338580706876, -0.183137063299809, -0.1143167151631 },
			std::initializer_list<double>{ -0.0408762886483087, 0.169471530184917, 0.172239202367153 },
			std::initializer_list<double>{ 0.204783837080904, 0.0691467768954059, 0.137616092244876 },
			std::initializer_list<double>{ -0.220932977424746, 0.209976278318135, -0.106234053454003 },
			std::initializer_list<double>{ -0.0856459776517156, -0.252626742374123, -0.1631410718876 },
			std::initializer_list<double>{ -0.130288819580702, -0.018587418572944, -0.114646746812767 }
		});*/
}