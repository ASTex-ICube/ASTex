#ifndef __BENCHMARKER__H__
#define __BENCHMARKER__H__

#include "patchProcessor.h"
#include "utils.h"

/////////////////////

class ContentExchangeBenchmarker : public ASTex::Benchmarker
{
public:
	ContentExchangeBenchmarker();

	using ImageType = ASTex::ImageRGBu8;

	//step 1: requirements//

	void setInput(const ImageType &input) {m_input=input;}
	void setOutputDirectories(const std::string &OldDirectory,
							  const std::string &RandomDirectory,
							  const std::string &GetisGIDirectory,
							  const std::string &PCTSDirectory);

	//step 2: options//

	void setPCTSArgumentsFilePath(const std::string &path) {m_pctsPath = path;}

	void setGenerateOld(bool b) {m_generateOld = b;}
	void setGenerateRandom(bool b) {m_generateRandom = b;}
	void setGenerateGetisGI(bool b) {m_generateGetisGI = b;}
	void setGeneratePCTS(bool b) {m_generatePCTS = b;}

	void setOutputSize(size_t width, size_t height);
	void setNbOutputs(unsigned nbOutputs) {m_nbOutputs = nbOutputs;}
	void setNbContentsPerPatch(unsigned int nbContentsPerPatch) {m_nbContentsPerPatch=nbContentsPerPatch;}

	//step 3

	void generate();

private:

	void _generateParametersForPP();

	ImageType m_input;

	void generateOldMethod();
	void generateRandomMethod();
	void generateGetisGIMethod();
	void generatePCTSMethod();

	bool m_generateOld;
	bool m_generateRandom;
	bool m_generateGetisGI;
	bool m_generatePCTS;

	//result directories

	std::string m_pctsPath;

	std::string m_OldDirectory;
	std::string m_RandomDirectory;
	std::string m_GetisGIDirectory;
	std::string m_PCTSDirectory;

	unsigned int m_nbContentsPerPatch;
	size_t m_outputWidth, m_outputHeight;
	unsigned int m_nbOutputs;

	ASTex::ContentExchange::PatchProcessor<ImageType> *m_pp;
};

ContentExchangeBenchmarker::ContentExchangeBenchmarker() :
	Benchmarker(),
	m_generateOld(true),
	m_generateRandom(true),
	m_generateGetisGI(false),
	m_generatePCTS(false),
	m_pctsPath(),
	m_OldDirectory(),
	m_RandomDirectory(),
	m_GetisGIDirectory(),
	m_PCTSDirectory(),
	m_nbContentsPerPatch(3),
	m_outputWidth(512),
	m_outputHeight(512),
	m_nbOutputs(1),
	m_pp(nullptr)
{}

void ContentExchangeBenchmarker::setOutputDirectories(const std::string &OldDirectory,
													   const std::string &RandomDirectory,
													   const std::string &GetisGIDirectory,
													   const std::string &PCTSDirectory)
{
	m_OldDirectory = OldDirectory + "/";
	ASTex::create_directory(m_root + m_OldDirectory);
	m_RandomDirectory = RandomDirectory + "/";
	ASTex::create_directory(m_root + m_RandomDirectory);
	m_GetisGIDirectory = GetisGIDirectory + "/";
	ASTex::create_directory(m_root + m_GetisGIDirectory);
	m_PCTSDirectory = PCTSDirectory + "/";
	ASTex::create_directory(m_root + m_PCTSDirectory);
}

void ContentExchangeBenchmarker::setOutputSize(size_t width, size_t height)
{
	m_outputWidth = width;
	m_outputHeight = height;
}

void ContentExchangeBenchmarker::generateOldMethod()
{
	for(unsigned i=0; i<m_nbOutputs; ++i)
	{
		_generateParametersForPP();
		m_pp->setSeed(i+1);
		m_pp->fullProcess_oldMethod();
		ASTex::Mipmap<ImageType> outputTexture = m_pp->generate();
		outputTexture.texture().save(m_root + m_OldDirectory + "output_" + std::to_string(i) + ".png");
		delete m_pp;
	}
}

void ContentExchangeBenchmarker::generateRandomMethod()
{
	for(unsigned i=0; i<m_nbOutputs; ++i)
	{
		_generateParametersForPP();
		m_pp->setSeed(i+1);
		m_pp->patches_initRandom(16);
		m_pp->contents_initDefault();
		m_pp->contents_initRandom();
		ASTex::Mipmap<ImageType> outputTexture = m_pp->generate();
		outputTexture.texture().save(m_root + m_RandomDirectory + "output_" + std::to_string(i) + ".png");
		delete m_pp;
	}
}

void ContentExchangeBenchmarker::generateGetisGIMethod()
{
	_generateParametersForPP();
	//TODO
	m_pp->generate();
	delete m_pp;
}

void ContentExchangeBenchmarker::generatePCTSMethod()
{
	for(unsigned i=0; i<m_nbOutputs; ++i)
	{
		_generateParametersForPP();
		m_pp->setSeed(i+1);
		m_pp->patches_initRandom(16);
		m_pp->contents_initDefault();
		m_pp->contents_initRandom();
		m_pp->contents_enhancePCTS(m_pctsPath);
		ASTex::Mipmap<ImageType> outputTexture = m_pp->generate();
		std::cout << "Generated texture number " << i << " by PCTS method" << std::endl;
		outputTexture.texture().save(m_root + m_PCTSDirectory + "output_" + std::to_string(i) + ".png");
		delete m_pp;
	}
}

void ContentExchangeBenchmarker::generate()
{
	if(m_generateOld)
		generateOldMethod();
	if(m_generateRandom)
		generateRandomMethod();
	if(m_generateGetisGI)
		generateGetisGIMethod();
	if(m_generatePCTS)
		generatePCTSMethod();

}

void ContentExchangeBenchmarker::_generateParametersForPP()
{
	m_pp = new ASTex::ContentExchange::PatchProcessor<ImageType>;
	m_pp->setTexture(m_input);
	m_pp->setFilteringMode(NO_FILTER);
	m_pp->setNbContentsPerPatch(m_nbContentsPerPatch);
	m_pp->setOutputSize(m_outputWidth, m_outputHeight);
}


#endif
