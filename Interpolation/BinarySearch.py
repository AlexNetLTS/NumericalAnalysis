""" 二分查找 """


def sortedsearch(lst, x):
    if x < lst[0] or x > lst[-1]:
        raise ValueError('出错了')
    mid = (len(lst) - 1) // 2
    if lst[mid] <= x <= lst[mid + 1]:
        return mid
    else:
        if x > lst[mid]:
            return sortedsearch(lst[mid: ], x)
        else:
            return sortedsearch(lst[: mid], x)


if __name__ == '__main__':
    lst = [1, 2, 3]
    print(sortedsearch(lst, 2.5))
    