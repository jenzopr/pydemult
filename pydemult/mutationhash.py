import itertools

# With input from Valentine Svensson<valentine@nxn.se

def mutationhash(strings, nedit, alphabet = ["A", "C", "G", "T"], log = None):
    """
    produce a hash with each key a nedit distance substitution for a set of
    strings. values of the hash is the set of strings the substitution could
    have come from
    """
    maxlen = max([len(string) for string in strings])
    indexes = generate_idx(maxlen, nedit, alphabet = alphabet)
    muthash = dict()
    for string in strings:
        if string not in muthash:
            if log is not None:
                log.debug("Added {} -> {} to mutationhash".format(string, string))
            muthash[string] = {string}
        for x in substitution_set(string, indexes):
            if x not in muthash:
                muthash[x] = {string}
                if log is not None:
                    log.debug("Added {} -> {} to mutationhash ".format(x, string))
            else:
                muthash[x].add(string)
                if log is not None:
                    log.debug("Added {} -> {} to mutationhash ".format(x, string))
    return muthash

def substitution_set(string, indexes):
    """
    for a string, return a set of all possible substitutions
    """
    strlen = len(string)
    return [mutate_string(string, x) for x in indexes if valid_substitution(strlen, x)]

def valid_substitution(strlen, index):
    """
    skip performing substitutions that are outside the bounds of the string
    """
    values = index[0]
    return all([strlen > i for i in values])

def generate_idx(maxlen, nedit, alphabet = ["A", "C", "G", "T"]):
    """
    generate all possible nedit edits of a string. each item has the form
    ((index1, index2), 'A', 'G')  for nedit=2
    index1 will be replaced by 'A', index2 by 'G'

    this covers all edits < nedit as well since some of the specified
    substitutions will not change the base
    """
    indexlists = []
    ALPHABETS = [alphabet for x in range(nedit)]
    return list(itertools.product(itertools.combinations(range(maxlen), nedit),
                                  *ALPHABETS))

def mutate_string(string, tomutate):
    strlist = list(string)
    for i, idx in enumerate(tomutate[0]):
        strlist[idx] = tomutate[i+1]
    return "".join(strlist)
