from functools import lru_cache
import time

@lru_cache(None)
def f(n):
    if n < 9: return n // 3 + n % 3
    if n >= 9: return f(n // 9) + f(n % 9)

count = 0
for i, n in enumerate(range(9**9, 0, -1)):
    # t = time.time()
    # if i<1000:
    #     print(n, f(n))
    if f(n)==33:
        count += 1
    # print(time.time()-t)
print(count)