#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''

This script tokenizes documents indexed in a  
file, normalizes tokens, creates a LDA-C formated 
corpus  

Created by: Clint P. George
Created On: Dec 14, 2013   

'''
import os 
import re 
import codecs
import logging

from gensim import corpora
from nltk.stem.porter import PorterStemmer
from nltk.tokenize import PunktWordTokenizer, RegexpTokenizer
from nltk.stem.wordnet import WordNetLemmatizer
    
    
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

'''
Global variables 
'''
STRIP_CHAR_LIST = [u'_', u'-', u',', u'!', u':', u'.', u'?', 
                   u';', u'=', u'…', u'•', u'–', u'¿', u'¡', 
                   u'º', u'ª', u'«', u'»', u'*', u'~', u'`', 
                   u':', u'<', u'>', u'{', u'}',
                   u'[', u']', u'//', u'(', u')'] 
REMOVE_LIST = [u"[", u"]", u"{", u"}", u"(", u")", u"'", u".", 
               u"..", u"...", u",", u"?", u"!", u"/", u"\"", 
               u"\"", u";", u":", u"-", u"`", u"~", u"@", u"$", 
               u"^", u"|", u"#", u"=", u"*"];

numeric_parser = re.compile(r"""        # A numeric string consists of:
    (?P<sign>[-+])?              # an optional sign, followed by either...
    (
        (?=\d|\.\d)              # ...a number (with at least one digit)
        (?P<int>\d*)             # having a (possibly empty) integer part
        (\.(?P<frac>\d*))?       # followed by an optional fractional part
        (E(?P<exp>[-+]?\d+))?    # followed by an optional exponent, or...
    |
        Inf(inity)?              # ...an infinity, or...
    |
        (?P<signal>s)?           # ...an (optionally signaling)
        NaN                      # NaN
        (?P<diag>\d*)            # with (possibly empty) diagnostic info.
    )
    \Z
""", re.VERBOSE | re.IGNORECASE | re.UNICODE).match
    
date_parser = re.compile(r'''(\d{2})[/.-](\d{2})[/.-](\d{4})$''',
                         re.VERBOSE | re.IGNORECASE | re.UNICODE).match
time_parser = re.compile(r'''(\d{2}):(\d{2}):(\d{2})$''',
                         re.VERBOSE | re.IGNORECASE | re.UNICODE).match
time_parser2 = re.compile(r'''(\d{2}):(\d{2})$''',
                         re.VERBOSE | re.IGNORECASE | re.UNICODE).match
'''
Initializes the tokenizer, lemmatizer, and stemmer  
'''
# Manages real numbers, numbers with comma separators, 
# dates, times, hyphenated words, email addresses, urls 
# and abbreviations
# https://code.google.com/p/nltk/issues/detail?id=128
pat1 = r'''(?:[A-Z][.])+|\d[\d,.:\-/\d]*\d|\w+[\w\-\'.&|@:/]*\w+'''

pat2 = r'''\w+|\$[\d\.]+|\S+'''

# https://code.google.com/p/nltk/issues/detail?id=128
pat3 = r'''(?x)    # set flag to allow verbose regexps
        ([A-Z]\.)+        # abbreviations, e.g. U.S.A. 
        | \w+(-\w+)*        # words with optional internal hyphens
        | \$?\d+(\.\d+)?%?  # currency and percentages, e.g. $12.40, 82%
        | \.\.\.            # ellipsis
        | [][.,;"'?():-_`]  # these are separate tokens
        '''

regx_tokenizer = RegexpTokenizer(pat1)

punkt_tokenizer = PunktWordTokenizer()
wordnet_lmtzr = WordNetLemmatizer()
porter_stemmer = PorterStemmer() 

def xstr(s):
    return '' if s is None else str(s)


def load_en_stopwords(filename):
    '''Loads English stop-words from a given file 
    
    Return: 
        a list of stop words
    Arguments: 
        the stop-words file name
    '''
    
    stopwords = []
    with codecs.open(filename, mode='r', encoding='utf-8') as fSW: 
        for line in fSW: 
            stopwords.append(line.strip().lower())
    return stopwords

stop_words = load_en_stopwords('en_stopwords')
print 'Number of stop words:', len(stop_words)

def cleanup(token):
    '''Clean up a given token based on regular expression, strip(),
    and a predefined set of characters 
    
    Returns: 
        a cleaned up token 
    Arguments:
        a token (str) 
    '''
    
    def remove_non_ASCII(s): 
        return u"".join(i for i in s if ord(i) < 128)
    
    try:
        
        token = remove_non_ASCII(token)

        for each_char in STRIP_CHAR_LIST:
            token = token.strip(each_char)
        
        if ((token in REMOVE_LIST) # removes if it's in the remove list 
            or (token in stop_words) # removes stop words  
            or numeric_parser(token) # discards if it's a numeric 
            or date_parser(token) # discards if it's a date  
            or time_parser(token) or time_parser2(token)): # discards if it's a time   
            return u""

    except: 
        return u"" 
    
    return token.strip()

def regex_tokenizer(text):
    '''A tokenizer based on NLTK's RegexpTokenizer 
    
    Returns: 
        a list of tokens 
    Arguments:
        a string to tokenized 
    
    '''

    tokens = regx_tokenizer.tokenize(text)

    filtered = []
    for w in tokens:
        try:
            token = cleanup(w.lower()) 
            if len(token) > 0: 
                filtered.append(token)
        except: 
            pass 

    return filtered

   
   
def punkt_word_tokenizer(text):
    '''A tokenizer based on NLTK's PunktWordTokenizer 
    
    Returns: 
        a list of tokens 
    Arguments:
        a string to tokenized 
    
    '''
    try: 
        text = ' '.join(text.lower().split()) # removes newline, tab, and white space
    except Exception:
        pass  
          
    tokens = punkt_tokenizer.tokenize(text)
    
    print tokens 
    # tokens = [cleanup(w) for w in tokens]
    # tokens = [w for w in tokens if w not in REMOVE_LIST]
    filtered = []
    for w in tokens:
        try: 
            w = cleanup(w)
            if len(w) > 0: 
                for wt in w.split(','): 
                    filtered.append(wt)
        except: 
            pass 
    
    return filtered

def whitespace_tokenize(doc_text):
    '''
    This function will tokenize the given 
    document text into tokens  based on whitespace 
    
    '''
    return doc_text.lower().split()



def lemmatize_tokens(word_tokens):
    '''
    Lemmatize tokens based on WordNet 
    '''
    tokens = [] 
    for token in word_tokens:
        try:
            # print 'lemma:', token, '-->', wordnet_lmtzr.lemmatize(token)   
            token = wordnet_lmtzr.lemmatize(token)
        except: pass 
        tokens.append(token)
    
    return tokens


def stem_tokens(word_tokens):
    '''
    Stem tokens based on Snowball stemmer  
    '''
    tokens = [] 
    for token in word_tokens:
        try:
            # print 'stem:', token, '-->', porter_stemmer.stem(token)   
            token = porter_stemmer.stem(token)
        except: pass 
        tokens.append(token)
        
    return tokens



                    


def build_lda_corpus(doc_path_index_file, data_folder,
                     dictionary_file, ldac_file, 
                     min_word_freq=5, 
                     min_word_len=2, max_word_len=20):
    '''
    The main function that does the job! 
    
    '''

    def process_index_doc(doc):
        '''Processes a single indexed document  
        
        Arguments: 
            doc - a document's details in the index file 
        '''
        
        _, _, _, doc_path = doc
        # print doc_id, category, subject, doc_path
        
        for body_charset in 'US-ASCII', 'ISO-8859-1', 'UTF-8':
            try:
                with codecs.open(doc_path, 'r', body_charset) as fp:
                    doc_text = fp.read()
                    doc_text = doc_text.encode('UTF-8') # encodes to UNICODE 
            except UnicodeError: 
                pass
            else: 
                break
        
        tokens = regex_tokenizer(doc_text)
        tokens = lemmatize_tokens(tokens)
        # tokens = stem_tokens(tokens)
        
        tokens = [w for w in tokens if w not in stop_words] # double checking 
            
        return tokens
    
    def create_dictionary(doc_details, 
                          dictionary_file, 
                          min_word_freq, 
                          min_word_len, max_word_len):
        '''
        Creates a dictionary from the given text files 
        using the Gensim class and functions
        
        Returns:
            None 
        Arguments:
            doc_details - a list of documents to be processed  
            dictionary_file - the dictionary object file (output)
            min_word_freq - min frequency of a valid vocabulary term 
            min_word_len - min word length of a valid vocabulary term 
            max_word_len - max word length of a valid vocabulary term
        '''
    
        # collect statistics about all tokens
        dictionary = corpora.Dictionary(process_index_doc(doc) 
                                            for doc in doc_details) # doc_id, category, subject, doc_path
        
        # remove words that appear < MIN_FREQUENCY 
        # remove words that have word length < MIN_WORD_LENGTH
        # remove words that have word length > MAX_WORD_LENGTH    
        filter_ids = [] 
        filter_ids += [tokenid for tokenid, docfreq in dictionary.dfs.iteritems() 
                       if docfreq < min_word_freq]
        filter_ids += [tokenid for tokenid, docfreq in dictionary.dfs.iteritems() 
                       if (len(dictionary[tokenid]) < min_word_len)]
        filter_ids += [tokenid for tokenid, docfreq in dictionary.dfs.iteritems() 
                       if (len(dictionary[tokenid]) > max_word_len)]
        
        dictionary.filter_tokens(filter_ids) # remove filtered words 
        dictionary.compactify() # remove gaps in id sequence after words that were removed
        dictionary.save(dictionary_file) # store the dictionary, for future reference
        
        logging.info(str(dictionary))
    
    
    class TextCorpus(object):
        '''The text corpus class 
        '''
        
        def __init__(self, dictionary, doc_details):
            '''
            Constructor 
            
            Returns: 
                a corpus object 
            Arguments: 
                dictionary - a dictionary object 
                doc_details - the document's path index file
            '''
            
            self._doc_details = doc_details               
            self._dictionary = dictionary
             
        def __iter__(self):
    
            for doc in self._doc_details:
                tokens = process_index_doc(doc)
                yield self._dictionary.doc2bow(tokens)
    

    
    # Reads the docs to be processed 
    doc_details = []
    with open(doc_path_index_file) as fp: 
        for lc, doc_line in enumerate(fp): 
            if lc == 0: continue
            doc_id, category, subject = doc_line.strip().split(';')
            doc_path = os.path.join(data_folder, subject)
            doc_details.append([doc_id, category, subject, doc_path])
    assert len(doc_details) > 0
    

    # Creates the dictionary 
    create_dictionary(doc_details, dictionary_file, 
                      min_word_freq, min_word_len, max_word_len)
    
    # Creates the corpus 
    dictionary = corpora.Dictionary().load(dictionary_file)       
    corpus = TextCorpus(dictionary, doc_details) # doesn't load the corpus into the memory!
    corpora.BleiCorpus.serialize(ldac_file, corpus, id2word=dictionary)
    
    
    index_file = os.path.splitext(ldac_file)[0] + '.index'
    with open(index_file, 'w') as fw:
        print >>fw, "docid;category;subject;uniquewordcount;doclength" 
        for doc_id, doc in enumerate(corpus): 
            doc_length = sum(count for _, count in doc)
            doc_id, category, subject, _ = doc_details[doc_id]
            print >>fw, ";".join([str(doc_id), category, subject, str(len(doc)), str(doc_length)]) 
        
    logging.info('The Blei corpus building is completed.')
    
    
    

data_folder = '/home/cgeorge/workspace/whales_tires/docs-cleaned'  
doc_path_index_file = '/home/cgeorge/workspace/whales_tires/whales-tires-cleaned.csv'
dictionary_file = '/home/cgeorge/workspace/whales_tires/whales-tires.dict' 
ldac_file = '/home/cgeorge/workspace/whales_tires/whales-tires.ldac'

build_lda_corpus(doc_path_index_file, data_folder, dictionary_file, ldac_file)

