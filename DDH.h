#pragma once

#include <iostream>
#include <openssl/sha.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include "NTL/ZZ.h"
#include "NTL/vector.h"

using namespace std;
using namespace NTL;

enum UpdateType
{
    MODIFY = 0,
    INSERT,
    DELETE
};

class DDHPIR{
public:
    ZZ N,p;
    Vec<ZZ> DB;
    Vec<ZZ> h;
    ZZ blocknum;
    ZZ blocksize;

    DDHPIR(long numbits){

        GenPrime(p, numbits + 2);
        power(N, 2, numbits);

//        DB.SetLength(to_long(N));
//        h.SetLength(to_long(N));

        SqrRoot(blocksize, N);

        blocknum = blocksize;

        if (blocknum * blocksize < N)
        {
            blocknum = blocknum + 1;
        }
        DB.SetLength(to_long(blocknum * blocksize));
        h.SetLength(to_long(blocknum * blocksize));

        for (int i = 0; i < N; i++)
        {
            DB[i] = RandomBnd(2);
            h[i] = i + 2;
        }

        for(int i = to_long(N); i < blocknum * blocksize; i++)
        {
            DB[i] = 1;
            h[i] = 1;
        }
    }
};

void offline(Vec<ZZ>& h, Vec<ZZ>& DB, Vec<ZZ>& digest, unsigned char hs[SHA256_DIGEST_LENGTH], ZZ blocknum, ZZ blocksize, ZZ N, ZZ p);

void query(Vec<ZZ>& h, ZZ i, ZZ p, Vec<ZZ>& q, Vec<ZZ>& params, ZZ blocknum, ZZ blocksize);

void Answer(Vec<ZZ>& DB, Vec<ZZ>& q, Vec<ZZ>& a, ZZ blocksize, ZZ blocknum, ZZ p, ZZ N);

void Recover(Vec<ZZ>& dc, Vec<ZZ>& h, Vec<ZZ>& a, unsigned char hs[SHA256_DIGEST_LENGTH], ZZ blocknum, ZZ l, ZZ i, ZZ r, ZZ p,
             ZZ t, ZZ *di);

void Update(Vec<ZZ>& DB, Vec<ZZ>& h, Vec<ZZ>& digest, Vec<ZZ>& DB_, Vec<ZZ>& h_, Vec<ZZ>& digest_, UpdateType updateType, Vec<ZZ>& array, long* index, ZZ blocksize, ZZ blocknum, ZZ blocknumIncremental, ZZ N, ZZ p, unsigned char hs_[SHA256_DIGEST_LENGTH]);
