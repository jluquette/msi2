/* #define DEBUG_SPUTNIK 1 */


/*
  find repeats in fasta format seq file 
  allows for indels, returns score.  

  beta version.  caveat emptor.  

  chrisa  29-Jul-94

  chris abajian
  University of Washington
  Dept. of Molecular Biotechnology  FJ-20
  Fluke Hall, Mason Road
  Seattle WA 98195
*/

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <stdlib.h>

/* trivial defs */
#ifndef True
#define True 1
#endif
#ifndef False
#define False 0
#endif

typedef int Boolean;

/* size of buffer for reads. */ 
#define BUF_SIZE 1024*10   /* 10K */
/* max size of description line (begins with ">") */
#define MAX_DESCRIPTION_LEN 1024
/* max sequence length */
#define MAX_SEQUENCE_LEN 1024*800 /* 800K */
/* max number of sequence chars dumped to line */
#define MAX_OUT_LINE_CHARS 60

/* for debugging only */
#define MAX_ERRCODES 1024

/* search params and definitions */
#define MIN_UNIT_LENGTH 1 /* start search with dinucleotide repeats */
/* will search for di, tri, tetra ... <n>nucleotide repeats up to 
   this value for n */
#define MAX_UNIT_LENGTH 4  /* up to and including pentanucleotides */
/* this is the point score for each exact match */
#define EXACT_MATCH_POINTS 1
/* this is the point score for a mismatch, insertion or deletion */
#define ERROR_MATCH_POINTS -6
/* this is the minimum score required to be considered a match */
#define MATCH_MIN_SCORE 8
/* this is the low score at which we stop trying */
#define MATCH_FAIL_SCORE -1
/* this is the max recursion depth we try to recover errors */
#define MAX_RECURSION 5


char *repeatName[MAX_UNIT_LENGTH+1] =
{
   "***ERROR***",    /* bad programmer!  no latte! */
   "mono",
   "di",
   "tri",
   "tetra"
};


char readBuf[BUF_SIZE];
Boolean endOfFile;
int curBufLen;
int curBufPos;
int fd;
Boolean havePutBack;
char putBack;

/* struct for indiv sequence in a file */
typedef struct ss
{
   char descStr[MAX_DESCRIPTION_LEN];
   char seqStr[MAX_SEQUENCE_LEN];
   unsigned int seqLen;
} SeqStruct, *SeqStructPtr;


/*
 * this structure describes the current state of a comparison.
 * it gets passed down to recursive calls of the find repeat
 * call so it can know when to bail out of an unsuccessful
 * search, or return the size/state of a successful hit, etc.
 */
typedef struct ms
{
   int curPos;         /* putative pattern starts here */
   int testPos;        /* start testing here */
   int testLen;        /* di, tri, tetra, etc. */
   int testCtr;        /* # chars in testLen already tested. mod counter */
   int curScore;       /* current score */
   int missense;       /* keep track of ins, del, err */
   int insertions;     
   int deletions;
   int depth;          /* how deep is recursion for this match */
   char errCodes[MAX_ERRCODES];
} MatchStruct, *MatchStructPtr;
/* a utility macro to copy one testStruct to another */
#define copyMSPtr(dest,source) memcpy((char *)dest,(char *)source,sizeof(MatchStruct))
/* a utility macro to increment the modular testCtr */
#define bumpTestCtr(msp) (msp)->testCtr++; if ((msp)->testCtr==(msp)->testLen) (msp)->testCtr=0;


/*
 ************************************************************
 * these routines are used to read and parse the fasta format
 * sequence file
 ************************************************************
 */

void fillBuf()
{
   size_t result;

   result = read(fd, (void *)readBuf, BUF_SIZE);
   if (result == -1)
     {
        fprintf(stderr,"error reading file! errno = %d\n",errno);
        exit(1);
     }
   else if (result == 0)
     endOfFile = True;
   else
     {
        curBufLen = result;
        curBufPos = 0;
     }
}  /* readBuf */


/* returns True on success */
Boolean getChar(char *achar)
{
   if (havePutBack)
     {
        *achar = putBack;
        havePutBack = False;
        return(True);
     }

   if (curBufPos == curBufLen)
     fillBuf();

   if (endOfFile)
     return (False);

   *achar = readBuf[curBufPos++];
   return (True);
}


void putCharBack(char c)
{
   havePutBack = True;
   putBack = c;
}


void openFile(char *fn)
{
   /* open the specified file */
   fd = open(fn, O_RDONLY);
   if (fd == -1)
     {
        fprintf(stderr,"unable to open file %s\n", fn);
        exit(1);
     }
}

/* should call this once for each file read */
void initBuffer()
{
   /* initialize length and pointer */
   curBufPos = 0;
   curBufLen = 0;
   havePutBack = False;
   endOfFile = False;
}

void addCharToLine(char c, char *line, int *lineLen)
{
   if (*lineLen < MAX_DESCRIPTION_LEN)
     line[(*lineLen)++] = c;
   else
     fprintf(stderr,"warning: description line truncated\n");
}


/*
 *********************************************************************
 * these routines are (more) specific to reading the fasta file format
 *********************************************************************
 */


/* 
 * pick up a non-blank line from the file, presumably description.
 * truncates all leading blanks and/or blank lines 
 */
Boolean getNonBlankLine(char *line)
{
   Boolean stop, nonBlank;
   char c;
   int lineLen;

   lineLen = 0;
   stop = False;
   nonBlank = False;  /* will be set by any non whitespace char */
   while ((! endOfFile) && (! stop))
     if (getChar(&c))
       if (c == '\n')
         stop = nonBlank; /* stop if have anything. don't save eol char. */
       else
         if (nonBlank)
           /* add it to line no matter what */
           addCharToLine(c,line,&lineLen);
         else if ((c != ' ') && (c != '\t'))
           {
              /* only non whitespace will start the line */
              nonBlank = True;
              addCharToLine(c,line,&lineLen);
           }

   return True;
}


/* load the sequence struct with comment line and bases */
SeqStructPtr getSeq(char *fname)
{
   SeqStructPtr newSeqP;
   Boolean endOfSeq;
   char c;

   if (endOfFile) return ((SeqStructPtr )0);   /* bombproofing */

   /* malloc a new seq */
   if (! (newSeqP = (SeqStructPtr )malloc(sizeof(SeqStruct)) ) )
     {
        fprintf(stderr,"unable to malloc() memory for sequence.\n");
        exit(1);
     }
   /* clear mem */
   memset( (void *)newSeqP, '\0', sizeof(SeqStruct));

   /* pick up description line */
   if (! getNonBlankLine(newSeqP->descStr) )
     {
        free(newSeqP);
        return ((SeqStructPtr )0);   
     }

   /* did it start correctly ? */
   if (newSeqP->descStr[0] != '>')
     {
        fprintf(stderr,"format error in input file:  missing '>'\n");
        exit(1);
     }

   endOfSeq = False;
   while ((!endOfFile) && (!endOfSeq))
     {
        if (getChar(&c))
          {
             if (c == '>')
               {
                  /* hit new sequence */
                  endOfSeq = True;
                  putCharBack(c);
               }
             else if (((c >= 'A') && (c <= 'Z')) ||
                      ((c >= 'a') && (c <= 'z')))/* bogus test, chris */
               /* have nucleotide */
               newSeqP->seqStr[newSeqP->seqLen++] = toupper(c);
             else if ((c != '\n') && (c != ' ') && (c != '\t'))
               {
                  /* wierd shit in file.  bail. */
                  fprintf(stderr,"bad char in sequence, file %s: %c\n",fname,c);
                  exit(1);
               } 
          }
     }     

   if (! newSeqP->seqLen)
     {
        fprintf(stderr,"? Null sequence encountered in file %s (ignored)\n",fname);
        fprintf(stderr,"  %s\n", newSeqP->descStr);
        free(newSeqP);
        return ((SeqStructPtr )0);   
     }

   return(newSeqP);
}  /* getSeq */


/* for debugging.  dump entire seq to stdout. */
#ifdef DEBUG_SPUTNIK
void dumpSeq(SeqStructPtr seqP)
{
   int i, charsOnLine;

   fprintf(stdout,"%s\n", seqP->descStr);
   fprintf(stdout,"Sequence (length = %d):\n", seqP->seqLen);
   i = 0;
   charsOnLine = 0;
   while (i < seqP->seqLen)
     {
        if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
        fprintf(stdout,"%c", seqP->seqStr[i++]);
     } 
   fprintf(stdout,"\n");
} /* dumpSeq */
#endif /* DEBUG_SPUTNIK */

/* dump the matched seq & stats to stdout */
void dumpMatch(SeqStructPtr seqP, 
               MatchStructPtr matchP,
               Boolean anyMatchThisSeq)
{
   int i, charsOnLine;

   //if (! anyMatchThisSeq)
     fprintf(stdout,"%s ", seqP->descStr);

   fprintf(stdout,"%s %d %d %d ",
           repeatName[matchP->testLen], 
           matchP->curPos+1,
           matchP->testPos,
           //matchP->testPos - matchP->curPos,
           matchP->curScore);

#ifdef DEBUG_SPUTNIK
   fprintf(stdout,"mis = %d, del = %d, ins = %d\n", 
           matchP->missense,
           matchP->deletions,
           matchP->insertions);
#endif

   //i = matchP->curPos;
   i = 0;
   charsOnLine = 0;
   //while (i < matchP->testPos)
   while (i < seqP->seqLen)
     {
        /*if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
		  */
        fprintf(stdout,"%c", seqP->seqStr[i++]);
     } 
   fprintf(stdout,"\n");

#ifdef DEBUG_SPUTNIK
   i = 0;
   charsOnLine = 0;
   while (i < (matchP->testPos - matchP->curPos))
     {
        if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
        if (matchP->errCodes[i] == '\0')
          fprintf(stdout," ");
        else
          fprintf(stdout,"%c", matchP->errCodes[i]);
        i++;
     } 
   fprintf(stdout,"\n");
#endif
}  /* dumpMatch */


Boolean testForNRepeat(SeqStructPtr seqP,
                       MatchStructPtr matchP)
{
   MatchStruct curMatch, recover, bestSoFar, bestOfABadLot;

   /* save matchP in case we fail altogether. */
   copyMSPtr(&curMatch, matchP);
   /* keep track of the best score and return that if over thresh. */
   copyMSPtr(&bestSoFar, matchP);

   while ( (curMatch.testPos < seqP->seqLen)           /* anything to test */
          && (curMatch.curScore > MATCH_FAIL_SCORE) )  /* above fail threshold */
     {
        /* test a base */
        if (seqP->seqStr[curMatch.curPos+curMatch.testCtr] 
            == seqP->seqStr[curMatch.testPos])
          {
             /* we matched.  this is easy. */
             curMatch.curScore += EXACT_MATCH_POINTS;  /* score your points */
             curMatch.testPos++; /* advance the downstream test position */
             bumpTestCtr(&curMatch); /* advance pos in the (presumed) repeating seq */
          }
        else if (seqP->seqStr[curMatch.testPos] == 'N')
          {
             /* don't call it wrong, but no credit either */
             curMatch.testPos++; /* advance the downstream test position */
             bumpTestCtr(&curMatch); /* advance pos in the (presumed) repeating seq */
          }
        else
          {
             /* no match.  take the score penalty, but keep going (maybe). */
             curMatch.curScore += ERROR_MATCH_POINTS;
             curMatch.testPos++; /* advance the downstream test position */
             bumpTestCtr(&curMatch);  /* advance pos in seq */
             /* is the score too bad to continue, or are we
                already too deep? */
             if ( (curMatch.curScore > MATCH_FAIL_SCORE)
                  && (curMatch.depth < MAX_RECURSION) )
               {
                  /* try simple missense */
                  copyMSPtr(&recover,&curMatch);
                  if ((recover.testPos - recover.curPos) < MAX_ERRCODES)
                    recover.errCodes[recover.testPos - recover.curPos -1] = 'M';
                  recover.missense++;
                  recover.depth++;
                  (void )testForNRepeat(seqP,&recover);
                  copyMSPtr(&bestOfABadLot,&recover);

                  /* try deletion */
                  copyMSPtr(&recover,&curMatch);
                  if ((recover.testPos - recover.curPos) < MAX_ERRCODES)
                    recover.errCodes[recover.testPos - recover.curPos -1] = 'D';
                  recover.testPos--; /* DON'T advance downstream */
                  recover.deletions++;
                  recover.depth++;
                  (void )testForNRepeat(seqP,&recover);
                  if (recover.curScore > bestOfABadLot.curScore)
                    copyMSPtr(&bestOfABadLot,&recover);

                  /* try insertion */
                  copyMSPtr(&recover,&curMatch);
                  if ((recover.testPos - recover.curPos) < MAX_ERRCODES)
                    recover.errCodes[recover.testPos - recover.curPos -1] = 'I';
                  /* RETEST for this base in the repeating seq */
                  if (recover.testCtr == 0)
                    recover.testCtr = recover.testLen - 1;
                  else
                    recover.testCtr--;
                  recover.insertions++;
                  recover.depth++;
                  (void )testForNRepeat(seqP,&recover);
                  if (recover.curScore > bestOfABadLot.curScore)
                    copyMSPtr(&bestOfABadLot,&recover);

                  /* take the best of a bad lot */
                  bestOfABadLot.depth--;  /* dec recursion counter */ 
                  copyMSPtr(&curMatch, &bestOfABadLot);
               }  /* it was worth carrying on */
          }  /* no match, found best of bad lot */

        /* whatever happened, the best we could do is now in matchP */
        if (curMatch.curScore > bestSoFar.curScore)
          copyMSPtr(&bestSoFar, &curMatch);

     }  /* while loop to test a single base */

   /* for whatever reason, we've stopped searching for more of this
      putative repeat.  if there were any matches that passed
      the global threshold, return the best of them. note that this
      has the effect of NOT advancing the pointer(s) if nothing
      rang the bell.  remember that we will test the same position
      for ntide repeats of several different lengths. */
   if (bestSoFar.curScore > MATCH_MIN_SCORE)
     {
        copyMSPtr(matchP, &bestSoFar);
        return(True);
     }
   return(False);   /* the whole thing was a waste of time */
}  /* testForNRepeat */


/* 
 * returns True if the sequence we want to look for repeats of is
 *
 * a) all the same base (i.e. 'AAA' or 'GG').  This filters out
 *    single nucleotide repeats
 *
 * b) conains 'N'.  we search against these, but don't use them
 *    as wildcards.
 */
Boolean ignoreSeq(SeqStructPtr seqP,
                  MatchStructPtr matchP)
{
   int i;

   /* firstly, never search for any pattern that contains N */
   for (i = 0; i < matchP->testLen; i++)
     if (seqP->seqStr[matchP->curPos+i] == 'N')
       return(True);

   /* now test for mononucleotide repeat.  other tests may get 
      added, in which case this one will beed to be changed. */
   // IGNORE this section to rescure MONONUCLEOTIDE repeats
   //for (i = 1; i < matchP->testLen; i++)
   //   if (seqP->seqStr[matchP->curPos] != seqP->seqStr[matchP->curPos+i])
   //     return(False);  /* they're not all the same */
   //return (True);  /* they ARE all same */
 
   return False;
}


void findRepeats(SeqStructPtr seqP)
{
   int curPos;
   Boolean anyMatchThisSeq, matchAtThisPos;
   MatchStruct match;

   memset( (char *)&match, 0, sizeof(MatchStruct) );  /* clear match struct */

   anyMatchThisSeq = False; /* avoid dumping description more than once. */
   /* loop on all positions in the sequence.  note that a match
      will advance curPos past all matching chars to the first
      unmatched char. */
   while ( match.curPos <= seqP->seqLen)
     {
        /* now loop on all the different lengths of repeats we're
           looking for (i.e. di, tri, tetra nucleotides.  if we
           find a match at a shorter repeat length, forego testing
           for longer lengths. */
        match.testLen = MIN_UNIT_LENGTH;
        matchAtThisPos = False;
        while ((match.testLen <= MAX_UNIT_LENGTH) && (!matchAtThisPos))
          {
             /* initialize the state of the match */
             match.curScore = 0;  /* no points yet */
             match.testCtr = 0; /* no chars tested yet */
             match.testPos = match.curPos + match.testLen;
             match.insertions = 0;
             match.deletions = 0;
             match.missense = 0;
             /* there are some things we don't want to test for (IGNORE - ignoring single repeat) */
             if (! ignoreSeq(seqP,&match))
               matchAtThisPos = testForNRepeat(seqP, &match);
             else
               matchAtThisPos = False;
             if (! matchAtThisPos) match.testLen++;
          }

        if (matchAtThisPos)
          {
             dumpMatch(seqP,&match,anyMatchThisSeq);
             anyMatchThisSeq |= matchAtThisPos;
             match.curPos = match.testPos;
          }
        else
          match.curPos++;  /* no, so advance to next base. */
     }
}


main(int argc, char* argv[])
{
   SeqStructPtr seqP;
   int count;

   if (argc != 2)
     {
        fprintf(stderr,"Usage: %s <fasta format sequence file name>\n", argv[0]);
        exit(1);
     }

   openFile(argv[1]);

   initBuffer();

   count = 0;
   while (! endOfFile) 
     if (seqP = getSeq(argv[1]))
       {
#ifdef DEBUG_SPUTNIK
          fprintf(stdout,"processing sequence %d\n", count++);
#endif
          /* dumpSeq(seqP); */
          findRepeats(seqP);
          free((void *)seqP);
       }
}
