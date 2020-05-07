#pragma once
#ifndef System_namespace
#define System_namespace

namespace System
{
	int I_BLOCKSIZE = 16;
	int J_BLOCKSIZE = 16;
	int K_BLOCKSIZE = 16;

	int count_it = 10;

	int getI_BLOCKSIZE() { return I_BLOCKSIZE; }
	int getJ_BLOCKSIZE() { return J_BLOCKSIZE; }
	int getK_BLOCKSIZE() { return K_BLOCKSIZE; }


	void static setJ_BLOCKSIZE(int blockSize) { J_BLOCKSIZE = blockSize; }
	void static setK_BLOCKSIZE(int blockSize) { K_BLOCKSIZE = blockSize; }
	void static setI_BLOCKSIZE(int blockSize) { I_BLOCKSIZE = blockSize; }
}
#endif //!end namespace