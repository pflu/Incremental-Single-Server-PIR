#include "DDH.h"

using namespace std;
using namespace NTL;

vector<unsigned char> ZZToByteArray(const ZZ& z)
{
    vector<unsigned char> bytes;
    ZZ temp = z;

    // 处理负数
    bool isNegative = (temp < 0);
    if (isNegative)
    {
        temp = -temp;
    }

    // 将 ZZ 类型的大整数转换为字节数组
    while (temp > 0)
    {
        unsigned long byte = conv<unsigned long>(temp % 256);
        bytes.push_back(static_cast<unsigned char>(byte));
        temp /= 256;
    }

    // 如果是负数，处理最高位的符号位
    if (isNegative)
    {
        for (size_t i = 0; i < bytes.size(); ++i)
        {
            bytes[i] = ~bytes[i];
        }
        // 加 1 以完成二进制补码表示
        for (size_t i = 0; i < bytes.size(); ++i)
        {
            if (++bytes[i] != 0)
            {
                break;
            }
        }
    }

    // 反转字节数组，使其从高位到低位排列
    reverse(bytes.begin(), bytes.end());

    return bytes;
}

void vectorcopy(ZZ* out, ZZ* data, long size){
    for(int i = 0; i < size; i++){
      out[i] = data[i];
    }
}

void DBUpdate(Vec<ZZ>& DB, Vec<ZZ>& out, UpdateType updateType, Vec<ZZ>& array, long* index, ZZ blocksize, ZZ blocknum, ZZ blocknumIncremental)
{
    long m = to_long(blocknumIncremental);
    if(updateType == MODIFY)
    {
        vectorcopy(out.elts(), DB.elts(), to_long(index[0] * blocksize));

        for(int i = 0; i < m - 1; i++)
        {
            vectorcopy(out.elts() + to_long(index[i] * blocksize), array.elts() + to_long(i * blocksize), to_long(blocksize));
            vectorcopy(out.elts() + to_long((index[i] + 1) * blocksize), DB.elts() + to_long((index[i] + 1) * blocksize), to_long((index[i + 1] - index[i] - 1) * blocksize));

        }
        vectorcopy(out.elts() + to_long(index[m - 1] * blocksize), array.elts() + to_long((m - 1) * blocksize), to_long(blocksize));
        vectorcopy(out.elts() + to_long((index[m - 1] + 1) * blocksize), DB.elts() + to_long((index[m - 1] + 1) * blocksize), to_long((blocknum - index[m - 1] - 1) * blocksize));
    }
    else if(updateType == INSERT)
    {
          vectorcopy(out.elts(), DB.elts(), to_long(blocksize * blocknum));
          vectorcopy(out.elts() + to_long(blocksize * blocknum), array.elts(), to_long(m * blocksize));
    }
    else if(updateType == DELETE)
    {
        vectorcopy(out.elts(), DB.elts(), to_long(index[0] * blocksize));
        for(int i = 0; i < m - 1; i++){
            vectorcopy(out.elts() + to_long((index[i] - i) * blocksize), DB.elts() + to_long((index[i] + 1) * blocksize), to_long((index[i + 1] - index[i] - 1) * blocksize));
        }
        vectorcopy(out.elts() + to_long((index[m - 1] - m + 1) * blocksize), DB.elts() + to_long((index[m - 1] + 1)*blocksize),to_long((blocknum - index[m - 1] - 1) * blocksize));
    }
}

void StateUpdate(Vec<ZZ>& h, Vec<ZZ>& h_, Vec<ZZ>& digest, Vec<ZZ>& digest_, Vec<ZZ>& D, UpdateType updateType, Vec<ZZ>& array, long* index, ZZ blocksize, ZZ blocknum, ZZ blocknumIncremental, ZZ N, ZZ p)
{
    if(updateType == MODIFY)
    {
        h_.SetLength(to_long(N));
        vectorcopy(digest_.elts(), digest.elts(), to_long(index[0]));
        for(int i = 0; i < blocknumIncremental - 1; i++){
            digest_[to_long(index[i])] = 1;
            for (int j = 0; j < blocksize; j++)
            {
                ZZ temp;

                PowerMod(temp, h[to_long(index[i] * blocksize + j)], D[to_long(index[i] * blocksize + j)], p);

                digest_[to_long(index[i])] = MulMod(digest_[to_long(index[i])], temp, p);
            }
            vectorcopy(digest_.elts() + index[i] + 1, digest.elts() + index[i] + 1, index[i + 1] - index[i] - 1);
        }
        vectorcopy(digest_.elts() + to_long(index[to_long(blocknumIncremental - 1)]), digest.elts() + to_long(index[to_long(blocknumIncremental - 1)]), to_long(blocknum - (index[to_long(blocknumIncremental - 1)] - 1)));
        vectorcopy(h_.elts(), h.elts(), to_long(blocksize * blocknum));
    }
    else if(updateType == INSERT)
    {
        h_.SetLength(to_long((blocknum + blocknumIncremental) * blocksize));
        vectorcopy(h_.elts(), h.elts(), to_long(blocksize * blocknum));
        for(int i = 0; i < blocknumIncremental * blocksize; i++)
        {
          h_[to_long(blocknum * blocksize) + i] = N + i + 2;
        }
        vectorcopy(digest_.elts(), digest.elts(), to_long(blocknum));
        for(int i = 0; i < blocknumIncremental; i++){
            digest_[to_long(blocknum + i)] = 1;
            for (int j = 0; j < blocksize; j++)
            {
                ZZ temp;

                PowerMod(temp, h_[to_long(blocknum * blocksize + i * blocksize + j)], D[to_long(blocknum * blocksize + i * blocksize + j)], p);

                digest_[to_long(blocknum + i)] = MulMod(digest_[to_long(blocknum + i)], temp, p);
            }
        }
    }
    else if(updateType == DELETE)
    {
        vectorcopy(digest_.elts(), digest.elts(), to_long(index[0]));
        vectorcopy(h_.elts(), h.elts(), to_long(index[0] * blocksize));
        for(int i = 0; i < blocknumIncremental - 1; i++)
        {
            vectorcopy(digest_.elts() + (index[i] - i), digest.elts() + (index[i] + 1), to_long(index[i + 1] - index[i] - 1));
            vectorcopy(h_.elts() + to_long((index[i] - i) * blocksize), h.elts() + to_long((index[i] + 1) * blocksize), to_long((index[i + 1] - index[i] - 1) * blocksize));
        }
        vectorcopy(digest_.elts() + to_long(index[to_long(blocknumIncremental - 1)] - blocknumIncremental + 1), digest.elts() + (index[to_long(blocknumIncremental - 1)] + 1), to_long(blocknum - index[to_long(blocknumIncremental - 1)] - 1));
        vectorcopy(h_.elts() + to_long((index[to_long(blocknumIncremental - 1)] - blocknumIncremental + 1) * blocksize), h.elts() + to_long((index[to_long(blocknumIncremental - 1)] + 1) * blocksize), to_long((blocknum - index[to_long(blocknumIncremental - 1)] - 1) * blocksize));
    }
}

void hash_SHA256(Vec<ZZ>& digest_, ZZ blocknum, ZZ blocknumIncremental, UpdateType updateType, unsigned char hs_[SHA256_DIGEST_LENGTH])
{
    if(updateType == MODIFY)
    {
        SHA256_CTX sha256;
        SHA256_Init(&sha256);
        for (int i = 0; i < blocknum; i++)
        {
            vector<unsigned char> temp = ZZToByteArray(digest_[i]);
            SHA256_Update(&sha256, temp.data(), temp.size());
        }

        SHA256_Final(hs_, &sha256);
    }
    if(updateType == INSERT)
    {
        SHA256_CTX sha256;
        SHA256_Init(&sha256);
        for (int i = 0; i < blocknum + blocknumIncremental; i++)
        {
            vector<unsigned char> temp = ZZToByteArray(digest_[i]);
            SHA256_Update(&sha256, temp.data(), temp.size());
        }

        SHA256_Final(hs_, &sha256);
    }
    if(updateType == DELETE)
    {
        SHA256_CTX sha256;
        SHA256_Init(&sha256);
        for (int i = 0; i < blocknum - blocknumIncremental; i++)
        {
            vector<unsigned char> temp = ZZToByteArray(digest_[i]);
            SHA256_Update(&sha256, temp.data(), temp.size());
        }

        SHA256_Final(hs_, &sha256);
    }
}

void Update(Vec<ZZ>& DB, Vec<ZZ>& h, Vec<ZZ>& digest, Vec<ZZ>& DB_, Vec<ZZ>& h_, Vec<ZZ>& digest_, UpdateType updateType, Vec<ZZ>& array, long* index, ZZ blocksize, ZZ blocknum, ZZ blocknumIncremental, ZZ N, ZZ p, unsigned char hs_[SHA256_DIGEST_LENGTH])
{
    DBUpdate(DB, DB_, updateType, array, index, blocksize, blocknum, blocknumIncremental);

    StateUpdate(h, h_, digest, digest_, DB_, updateType, array, index, blocksize, blocknum, blocknumIncremental, N, p);

    hash_SHA256(digest_, blocknum, blocknumIncremental, updateType, hs_);
}

void offline(Vec<ZZ>& h, Vec<ZZ>& DB, Vec<ZZ>& digest, unsigned char hs[SHA256_DIGEST_LENGTH], ZZ blocknum, ZZ blocksize, ZZ N, ZZ p)
{
    for (int i = 0; i < blocknum; i++)
    {
        digest[i] = 1;
        for (int j = 0; j < blocksize; j++)
        {
            ZZ temp;

            PowerMod(temp, h[to_long(i * blocksize + j)], DB[to_long(i * blocksize + j)], p);

           digest[i] = MulMod(digest[i], temp, p);
        }
    }
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    for (int i = 0; i < blocknum; i++)
    {
        vector<unsigned char> temp = ZZToByteArray(digest[i]);
        SHA256_Update(&sha256, temp.data(), temp.size());
    }

    SHA256_Final(hs, &sha256);
}

void query(Vec<ZZ>& h, ZZ i, ZZ p, Vec<ZZ>& q, Vec<ZZ>& params, ZZ blocknum, ZZ blocksize)
{
    ZZ l;
    ZZ templ;

    rem(templ, i, blocksize);
    div(l, i, blocksize);

    ZZ j = i % blocksize;

    ZZ r, t;
    RandomBnd(r, p);
    RandomBnd(t, p);

    ZZ temp;
    mul(temp, l, blocksize);
    for (int i = 0; i < blocksize; i++)
    {
        q[i] = h[to_long(temp + i)];
        PowerMod(q[i], q[i], r, p);
    }
    ZZ multemp;
    PowerMod(multemp, h[to_long(temp + j)], t, p);
    q[to_long(j)] = MulMod(q[to_long(j)], multemp, p);;

    params[0] = i;
    params[1] = l;
    params[2] = r;
    params[3] = t;
}

void Answer(Vec<ZZ>& DB, Vec<ZZ>& q, Vec<ZZ>& a, ZZ blocksize, ZZ blocknum, ZZ p, ZZ N)
{
    for (int i = 0; i < blocknum; i++)
    {
        a[i] = 1;
        for (int j = 0; j < blocksize; j++)
        {
            ZZ temp;

            PowerMod(temp, q[j], DB[to_long(i * blocksize + j)], p);
            a[i] = MulMod(a[i], temp, p);
        }
    }
}

void Recover(Vec<ZZ>& dc, Vec<ZZ>& h, Vec<ZZ>& a, unsigned char hs[SHA256_DIGEST_LENGTH], ZZ blocknum, ZZ l, ZZ i, ZZ r, ZZ p,
             ZZ t, ZZ *di)
{
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    for (int i = 0; i < blocknum; i++)
    {
        vector<unsigned char> temp = ZZToByteArray(dc[i]);
        SHA256_Update(&sha256, temp.data(), temp.size());
    }
    unsigned char hash[SHA256_DIGEST_LENGTH];

    SHA256_Final(hash, &sha256);

    for (int i = 0; i < SHA256_DIGEST_LENGTH; i++)
    {
      if (hs[i] != hash[i]){
        cout<<"Digest does not match"<<endl;
        //assert(false);
      }
    }

    ZZ tempdi = PowerMod(dc[to_long(l)], r, p);
    ZZ tempdiinv = InvMod(tempdi, p);
    tempdi = MulMod(tempdiinv, a[to_long(l)], p);

    if (tempdi == 1)
    {
        *di = 0;
    }
    else if (tempdi == PowerMod(h[to_long(i)], t, p))
    {
        *di = 1;
    }
    else
    {
      cout << "DDH PIR is false" << endl;
      assert(false);
    }
}