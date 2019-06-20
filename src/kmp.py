maxn = 1e5 + 10
def cal_next(s, plen):
    nxt = [0]*(max(plen,len(s))+2)
    for i in range(1,plen):
        j = nxt[i]
        while j>0 and s[i] != s[j]:
            j = nxt[j]
        if (s[j]==s[i]) :
            nxt[i+1] = j+1
        else:
            nxt[i+1] = 0
    return nxt

def kmp(str_f, ptr):
    slen = len(str_f)
    plen = len(ptr)
    if plen == 0 or slen == 0:
        return 0
    nxt = cal_next(ptr, plen )

    j = 0
    for i in range(slen):
        while( j>0 and j <len(ptr) and ptr[j] != str_f[i]):
            j = nxt[j]
        if j<len(ptr):
            if ptr[j] == str_f[i]:
                j += 1
    return j

if __name__ == "__main__":
    str1 = "aatat"
    str2 = "tataa"
    match_len = kmp(str1,str2)
    print("{},{}".format(str1[-match_len:],match_len))