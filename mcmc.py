# https://math.uchicago.edu/~shmuel/Network-course-readings/MCMCRev.pdf
# http://probability.ca/jeff/ftpdir/decipherart.pdf

import datetime
import numpy as np
import string

from tqdm import tqdm

# Alphabet
# alphabet = list(string.printable)
alphabet = [' '] + list(string.ascii_lowercase)


def encrypt(msg):
    msg = msg.lower()
    # Encryption
    enc_order = np.arange(len(alphabet))
    np.random.shuffle(enc_order)
    enc_a = [alphabet[enc_order[alphabet.index(c)]] for c in alphabet]
    # Encrypt message
    enc_msg = ''.join([enc_a[alphabet.index(c)] for c in msg])
    # Decryption
    dec_order = [enc_a.index(c) for c in alphabet]
    dec_a = [alphabet[dec_order[alphabet.index(c)]] for c in alphabet]
    # Decrypt message
    dec_msg = ''.join([dec_a[alphabet.index(c)] for c in enc_msg])
    return enc_msg, enc_order, enc_a, dec_order, dec_a


def show_enc(enc_order):
    for i in range(len(alphabet)):
        print(''.join([alphabet[i], '->', alphabet[enc_order[i]]]))


def show_dec(dec_order):
    l = []
    for i in range(len(alphabet)):
        l.append([alphabet[dec_order[i]], '<-', alphabet[i]])
    # Sort
    l = sorted(l)
    for a in l:
        print(''.join(a))


def make_transition_matrix(text_file, truncate_to=None):
    # Read all characters in text file
    chars = []
    with open(text_file) as f:
        for line in f:
            for c in line.lower():
                if c in alphabet:
                    chars.append(c)
    # Truncate chars
    if truncate_to:
        chars = chars[:truncate_to]
    # Init M
    M = np.ones((len(alphabet), len(alphabet)))
    # Fill M
    for c in tqdm(range(len(chars)-1)):
        i, j = alphabet.index(chars[c]), alphabet.index(chars[c+1])
        M[i, j] += 1
    # Return
    M = M/np.reshape(np.sum(M, 1), (len(M), 1))
    return M, chars


def mcmc(enc_msg, text_file, truncate_to=None, iters=10000):
    M, chars = make_transition_matrix(text_file, truncate_to)
    # Initial estimate of dec_order
    est_dec_order = np.arange(len(alphabet))
    np.random.shuffle(est_dec_order)
    # Calculate plausibility
    pl = 0
    for c in range(len(enc_msg)-1):
        pl += np.log(M[est_dec_order[alphabet.index(enc_msg[c])], est_dec_order[alphabet.index(enc_msg[c+1])]])
    for i in range(iters):
        # Randomly swap two symbols
        choices = np.random.choice(len(alphabet), 2)
        new_dec_order = list(est_dec_order)
        new_dec_order[choices[0]], new_dec_order[choices[1]] = new_dec_order[choices[1]], new_dec_order[choices[0]]
        # Calculate plausibility
        new_pl = 0
        for c in range(len(enc_msg)-1):
            new_pl += np.log(M[new_dec_order[alphabet.index(enc_msg[c])], new_dec_order[alphabet.index(enc_msg[c+1])]])
        # Choose if greater
        if new_pl > pl:
            # Accept
            est_dec_order = new_dec_order
            pl = new_pl
        else:
            # Flip a coin
            coin_toss = np.random.uniform()
            # Choose to keep or not
            if coin_toss < np.exp(new_pl - pl):
                # Accept
                est_dec_order = new_dec_order
                pl = new_pl
        # Check decrypted message
        dec_a = [alphabet[est_dec_order[alphabet.index(c)]] for c in alphabet]
        dec_msg = ''.join([dec_a[alphabet.index(c)] for c in enc_msg])
        print('{0:%Y%m%d_%H%M%S}'.format(datetime.datetime.now()), i, dec_msg)


# Does not work for short or repetitive messages!
# msg = "i love you baby"
# msg = "i love you baby i miss you so much i love you i love you i will always love you"

# Works for longer messages!
# msg = "enter hamlet ham to be or not to be that is the question whether tis nobler in the mind to suffer the slings and arrows of outrageous fortune or to take arms against a sea of troubles and by opposing end"
msg = "the algorithm continues trying to improve the current f by making random transpositions the coin tosses allow it to go to less plausible fs and keep it from getting stuck in local maxima of course the space of fs is huge or so why should this Metropolis random walk succeed to investigate this marc tried the algorithm out on a problem to which he knew the answer figure shows a well known section of shakespeares hamlet"

enc_msg, enc_order, enc_a, dec_order, dec_a = encrypt(msg)

# text_file = 'war_and_peace.txt'
# M, chars = make_transition_matrix(text_file, truncate_to=None)

mcmc(enc_msg, text_file, iters=50000)
