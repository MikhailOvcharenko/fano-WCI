import math
import multiprocessing
import gc

# List of weak partitions of a number n into k parts

def weak_compositions(n, k):
    if (n < 0) or (k < 0):
        return;
    elif (k == 0):
        if (n == 0):
            yield [];
        return;
    elif (k == 1):
        yield [n];
        return;
    else:
        for i in range(0, n + 1):
            for comp in weak_compositions(n - i, k - 1):
                yield [i] + comp;

# Various restrictions on weights and degrees of a Fano WCI
# described in [arXiv:2006.05666, Section 2]

def is_potentially_fano(weight_tuple, degree_tuple):
    return (sum(weight_tuple) - weight_tuple.count(1) < \
    sum(degree_tuple) + len(degree_tuple));

def is_weakly_quasismooth(weight_tuple, degree_tuple):
    if (degree_tuple[-1] < 2 * weight_tuple[-1]):
        return False;

    for i in range(len(degree_tuple)):
        if not (degree_tuple[-(i + 1)] > weight_tuple[-(i + 1)]):
            return False;
        
    return True;

def is_weakly_regular(weight_tuple, degree_tuple):
    return (math.prod(degree_tuple).mod(math.prod(weight_tuple)) == 0);

def is_not_linear_cone(weight_tuple, degree_tuple):
    for i in weight_tuple:
        for j in degree_tuple:
            if (i == j) : return False;
    return True;

def is_regular_pair(weight_tuple, degree_tuple):
    N = len(weight_tuple);
    weight_values = [];
    for i in range(N) : weight_values.append(set());
        
    for I in Subsets(range(N)):
        L = [weight_tuple[i] for i in I];
        l = len(L); g = gcd(L);
        if (g < 2) : continue;
        weight_values[l - 1].add(g);
    
    c = len(degree_tuple);
    degree_values = [];
    for j in range(c) : degree_values.append(set());
        
    for J in Subsets(range(c)):
        L = [degree_tuple[j] for j in J];
        l = len(L); g = gcd(L);
        if (g < 2) : continue;
        degree_values[l - 1].add(g);
    
    for i in range(N):
        weight_index = -1;
        if (weight_values[- (i + 1)] != set()) :
            weight_index = N - i - 1;
            break;

    for j in range(c):
        degree_index = -1;
        if (degree_values[- (j + 1)] != set()) :
            degree_index = c - j - 1;
            break;

    if (weight_index > degree_index) : return False;
    m = degree_index + 1;
    
    for k in range(m):
        for i in weight_values[k]:
            check = 0;
            for j in degree_values[k]:
                if (j.mod(i) == 0) : check = 1;
            if not check : return False;
    
    return True;

def check_restrictions(weight_tuple, degree_tuple):
    degree_list = list(degree_tuple);
    if ((degree_list == sorted(degree_list)) and \
        (degree_list.count(2) == 0) and \
        is_not_linear_cone(weight_tuple, degree_tuple) and \
        is_potentially_fano(weight_tuple, degree_tuple) and \
        is_weakly_quasismooth(weight_tuple, degree_tuple) and \
        is_weakly_regular(weight_tuple, degree_tuple) and \
        is_regular_pair(weight_tuple, degree_tuple)):
        return True;
    else:
        return False;
    
# Compute the list of generating families of
# Fano CI of variance r and codimension c
    
def bruteforce_CI(r, c):
    if (r <= 0) : raise TypeError('Variance should be positive.');
    if (c <= 0) : raise TypeError('Codimension should be positive.');

    weight_tuple = tuple((r + 2*c + 1)*[1]);
    delta_tuples = list(weak_compositions(r - 1, c));
    output = [];
            
    for delta_tuple in delta_tuples:
        degree_vector = vector(ZZ, c);
        for j in range(c):
            degree_vector[j] = weight_tuple[j] + delta_tuple[j] + 1;
        degree_vector[-1] += 1;
        degree_tuple = tuple(degree_vector);

        if check_restrictions(weight_tuple, degree_tuple):
            output.append([weight_tuple, degree_tuple]);

    for I in output : print(I);

def bruteforce_WCI_big_O1_routine(weight_tuple):
    output = [];
    delta_tuples = list(weak_compositions(var - 1, codim));
    
    for delta_tuple in delta_tuples:
        degree_vector = vector(ZZ, codim);
        for j in range(codim):
            degree_vector[j] = weight_tuple[j] + delta_tuple[j] + 1;
        degree_vector[-1] += 1;
        degree_tuple = tuple(degree_vector);

        if check_restrictions(weight_tuple, degree_tuple):
            weight_tuple_extended = (var + codim + 1)*[1] + list(weight_tuple);
            output.append([tuple(weight_tuple_extended), degree_tuple]);

    return output;

# Compute the list of generating families of
# Fano WCI of variance r and codimension c
# assuming that dim(|O_X(1)|) >= dim(X)

def bruteforce_WCI_big_O1(r, c):
    if (r <= 0) : raise TypeError('Variance should be positive.');
    if (c <= 0) : raise TypeError('Codimension should be positive.');

    # Building the input list
    weight_range = list(range(1, r + 2));
    weight_tuples = set();
    for T in UnorderedTuples(weight_range, c):
        if (T != c * [1]):
            weight_tuples.add(tuple(T))

    global var; var = r;
    global codim; codim = c;
        
    # Parallel computing
    output = [];
    if __name__ == '__main__':
        gc.enable()
        p = multiprocessing.Pool(multiprocessing.cpu_count());        
        for result in p.imap_unordered(bruteforce_WCI_big_O1_routine, \
                                       weight_tuples, chunksize=1):
            output += result;
        p.close();
        p.join();

    for I in output : print(I);
    
# Compute the list of generating families of
# Fano WCI of variance r and codimension c
# assuming that dim(|O_X(1)|) < dim(X)
    
def bruteforce_WCI_small_O1_routine(weight_tuple):
    output = [];
    last_weights = weight_tuple[-codim:];
    delta = sum(weight_tuple[:var - 1])
    delta_tuples = list(weak_compositions(delta, codim));
        
    for delta_tuple in delta_tuples:
        degree_vector = vector(ZZ, codim);
        for j in range(codim):
            degree_vector[j] = last_weights[j] + delta_tuple[j] + 1;    
        degree_vector[-1] += weight_tuple[var - 1];
        degree_tuple = tuple(degree_vector);

        if (weight_tuple.count(1) >= var) : continue;

        if check_restrictions(weight_tuple, degree_tuple):
            weight_tuple_extended = (codim + 1)*[1] + list(weight_tuple);
            output.append([tuple(weight_tuple_extended), degree_tuple]);
            
    return output;
            
def bruteforce_WCI_small_O1(r, c):
    if (r <= 0) : raise TypeError('Variance should be positive.');
    if (c <= 0) : raise TypeError('Codimension should be positive.');
    
    # Building the input list
    weight_range = range(1, (c + 2*r) + 1);
    weight_tuples = set();
    for T in UnorderedTuples(weight_range, c + r):
        if (T != (c + r) * [1]):
            weight_tuples.add(tuple(T));

    global var; var = r;
    global codim; codim = c;

    # Parallel computing
    output = [];
    if __name__ == '__main__':
        gc.enable()
        p = multiprocessing.Pool(multiprocessing.cpu_count());        
        for result in p.imap_unordered(bruteforce_WCI_small_O1_routine, \
                                       weight_tuples, chunksize=1):
            output += result;
        p.close();
        p.join();
        
    for I in output : print(I);

# Compute the list of generating families of
# Fano CI or WCI of variance r and codimension <= r.
    
def bruteforce_WCI(var):
    if (var < 0) : raise TypeError('Variance should be non-negative.');

    if (var == 0):
        print([(1,), ()]);
        return;
    
    for codim in range(1, var + 1):
        bruteforce_CI(var, codim);
        bruteforce_WCI_big_O1(var, codim);
        bruteforce_WCI_small_O1(var, codim);
