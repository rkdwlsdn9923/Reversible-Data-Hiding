
#include "global.h"

#include "rewrite_coeff4x4.h"
#include "elements.h"
#include "mb_access.h"
#include "memalloc.h"
#include "rtp.h"
#include "vlc.h"

static int RBSPtoEBSP(byte *NaluBuffer, unsigned char *rbsp, int rbsp_size);

int init_Rewrite(int initok)
{
	rw_bs = (RW_Bitstream*)malloc(sizeof(RW_Bitstream));
	if ( rw_bs == NULL ) no_mem_exit("FinalizeSpareMBMap: rw_bs");
	rw_bs->streamBuffer = (byte*)malloc(MAXRTPPAYLOADLEN);
	if ( rw_bs->streamBuffer == NULL ) no_mem_exit("FinalizeSpareMBMap: rw_bs->streamBuffer");
	rw_bs->bits_to_go  = 8;
	rw_bs->byte_pos    = 0;
	rw_bs->byte_buf    = 0;
	memset( rw_bs->streamBuffer, 0, MAXRTPPAYLOADLEN);
	initok=1;

	return initok;

}

void free_Rewrite(RW_Bitstream *rw_bs)
{
	if (rw_bs != NULL)
	{
		if (rw_bs->streamBuffer != NULL)
		{
			free(rw_bs->streamBuffer);       
			rw_bs->streamBuffer = NULL;
		}
		free(rw_bs);
		rw_bs = NULL;
	}
}

//buffer file에 쓰기(frame 별로)
int FlushStreamBuffer(RW_Bitstream *rw_bs, FILE *f_annexb)
{
	int Bufflen = rw_bs->byte_pos + (int)((double)rw_bs->byte_pos / (double)3);
	BYTE* Buffer;

	if ((Buffer = (byte*)calloc (Bufflen, sizeof (byte))) == NULL)
  {
    free (Buffer);
    no_mem_exit ("AllocBuffer: Buffer");
  }

	Bufflen = RBSPtoEBSP(Buffer, rw_bs->streamBuffer, rw_bs->byte_pos);

	if (Bufflen != fwrite (Buffer, 1,Bufflen, f_annexb))
  {
    printf ("Fatal: cannot write %d bytes to bitstream file, exit (-1)\n", Bufflen);
    exit (-1);
  }

	free(Buffer);
  return Bufflen;
}

/*!
************************************************************************
*  \brief
*     This function add emulation_prevention_three_byte for all occurrences
*     of the following byte sequences in the stream
*       0x000000  -> 0x00000300
*       0x000001  -> 0x00000301
*       0x000002  -> 0x00000302
*       0x000003  -> 0x00000303
*
*  \param NaluBuffer
*            pointer to target buffer
*  \param rbsp
*            pointer to source buffer
*  \param rbsp_size
*           Size of source
*  \return
*           Size target buffer after emulation prevention.
*
************************************************************************
*/

static int RBSPtoEBSP(byte *NaluBuffer, unsigned char *rbsp, int rbsp_size)
{
  int j     = 0;
  int count = 0;
  int i;

  for(i = 0; i < rbsp_size; i++)
  {
    if(count == ZEROBYTES_SHORTSTARTCODE && !(rbsp[i] & 0xFC))
    {
      NaluBuffer[j] = 0x03;
      j++;
      count = 0;
    }
    NaluBuffer[j] = rbsp[i];
    if(rbsp[i] == 0x00)
      count++;
    else
      count = 0;
    j++;
  }

  return j;
}



void RewirteCopyLevRun(int* levarr, int* runarr, int *rw_levarr,int *rw_runarr, int numcoef)
{
	memset(rw_levarr, 0, sizeof(int)*17);
	memset(rw_runarr, 0, sizeof(int)*17);

	memcpy(rw_levarr, levarr, sizeof(int)*numcoef);
	memcpy(rw_runarr, runarr, sizeof(int)*numcoef);
}


/*!
 ************************************************************************
 * \brief
 *    Writes coeff of an 4x4 block (CAVLC)
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */
int RewriteCoeff4x4_CAVLC_normal (Macroblock* currMB, int block_type, int b8, int b4, int *levarr, int *runarr, RW_Bitstream *bsBuff)
{
  Slice* currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currSlice->p_Vid;
  int           no_bits    = 0;
  SyntaxElement se;

  //RW_Bitstream *bs;
  //DataPartition *dataPart;
  const byte           *partMap   = assignSE2partition[currSlice->dp_mode];

  int k,level = 1,run = 0, vlcnum;
  int numcoeff = 0, lastcoeff = 0, numtrailingones = 0; 
  int numones = 0, totzeros = 0, zerosleft, numcoef;
  int numcoeff_vlc;
  int code, level_two_or_higher;
  int dptype = 0;
  int nnz, max_coeff_num = 0, cdc = 0, cac = 0;
  int subblock_x, subblock_y;
  //int *mb_bits_coeff = &currMB->bits.mb_y_coeff;
#if TRACE
  char type[15];
#endif

  static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768};  // maximum vlc = 6


  int*  pLevel = NULL;
  int*  pRun = NULL;

    static const unsigned char chroma_ac_param[3][8][4] =
  {
    {{ 4, 20,  5, 21},
    {36, 52, 37, 53},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0}},
    {{ 4, 20,  5, 21},
    { 6, 22,  7, 23},
    {36, 52, 37, 53},
    {38, 54, 39, 55},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0},
    { 0,  0,  0,  0}},
    {{ 4, 20,  5, 21},
    {36, 52, 37, 53},
    { 6, 22,  7, 23},
    {38, 54, 39, 55},
    { 8, 24,  9, 25},
    {40, 56, 41, 57},
    {10, 26, 11, 27},
    {42, 58, 43, 59}}
  };

  int   yuv = p_Vid->yuv_format - 1;

  switch (block_type)
  {
  case LUMA:
    max_coeff_num = 16;

    pLevel = levarr;
    pRun   = runarr;
	/*pLevel = currSlice->cofAC[b8][b4][0];
    pRun   = currSlice->cofAC[b8][b4][1];*/
#if TRACE
    sprintf(type, "%s", "Luma");
#endif
    dptype = (currMB->is_intra_block == TRUE) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
	//dptype = (is_intra (currMB)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    break;

  case CHROMA_AC:
    max_coeff_num = 15;
    //mb_bits_coeff = &currMB->bits.mb_uv_coeff;
    cac = 1;

    pLevel = levarr;
    pRun   = runarr;
	/*pLevel = currSlice->cofAC[b8][b4][0];
    pRun   = currSlice->cofAC[b8][b4][1];*/
#if TRACE
    sprintf(type, "%s", "ChrAC");
#endif
    dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
	//dptype = (is_intra (currMB)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    break;

  case CHROMA_DC:
    max_coeff_num = p_Vid->num_cdc_coeff;
    //mb_bits_coeff = &currMB->bits.mb_uv_coeff;
    cdc = 1;
	
	pLevel = levarr;
    pRun   = runarr;
	/*pLevel = currSlice->cofDC[param + 1][0];
    pRun   = currSlice->cofDC[param + 1][1];*/
#if TRACE
    sprintf(type, "%s", "ChrDC");
#endif
    dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
	//dptype = (is_intra (currMB)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    break;

  case LUMA_INTRA16x16AC:
    max_coeff_num = 15;

    pLevel = levarr;
    pRun   = runarr;
	/*pLevel = currSlice->cofAC[b8][b4][0];
    pRun   = currSlice->cofAC[b8][b4][1];*/
#if TRACE
    sprintf(type, "%s", "Lum16AC");
#endif
    dptype = SE_LUM_AC_INTRA;
    break;

  case LUMA_INTRA16x16DC:
    max_coeff_num = 16;

	pLevel = levarr;
    pRun   = runarr;
    /*pLevel = currSlice->cofDC[0][0];
    pRun   = currSlice->cofDC[0][1];*/
#if TRACE
    sprintf(type, "%s", "Lum16DC");
#endif 
    dptype = SE_LUM_DC_INTRA;

	//fastmemset(bsBuff, 0, sizeof());
    break;


  default:
    error("writeCoeff4x4_CAVLC: Invalid block type", 600);
    break;
  }
 
  //bs = rw_bs;
  //dataPart = &(currSlice->partArr[partMap[dptype]]);

  
  for(k = 0; (k <= ((cdc) ? p_Vid->num_cdc_coeff : 16)) && level != 0; k++)
  {
    level = pLevel[k]; // level
    run   = pRun[k];   // run

    if (level)
    {

      totzeros += run; 
      if (iabs(level) == 1)
      {
        numones ++;
        numtrailingones ++;
        numtrailingones = imin(numtrailingones, 3); // clip to 3
      }
      else
      {
        numtrailingones = 0;
      }
      numcoeff ++;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      lastcoeff = k;
    }
  }

  if (!cdc)
  {
    if (!cac)
    {
      // luma
      subblock_x = b8;
     // subblock_x = ((b8 & 0x1) == 0) ? (((b4 & 0x1) == 0) ? 0 : 1) : (((b4 & 0x1) == 0) ? 2 : 3);
      // horiz. position for coeff_count context
      subblock_y = b4;
      //subblock_y = (b8 < 2) ? ((b4 < 2) ? 0 : 1) : ((b4 < 2) ? 2 : 3);
      // vert.  position for coeff_count context
      nnz = Repredict_nnz(currMB, LUMA, subblock_x,subblock_y);
    }
    else
    {
      // chroma AC
		int param = chroma_ac_param[yuv][b8][b4];
		/*subblock_x = b8>>2;
		subblock_y = (b4-4);*/
      subblock_x = param >> 4;
      subblock_y = param & 15;
      nnz = Repredict_nnz_chroma(currMB, subblock_x, (subblock_y-4)*4);
    }
//    p_Vid->nz_coeff [currMB->mbAddrX ][0][subblock_y][subblock_x] = numcoeff;
	rdh_nz_coeff  [currMB->mbAddrX ][subblock_y][subblock_x] = numcoeff;

    numcoeff_vlc = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));
  }
  else
  {
    // chroma DC (has its own VLC)
    // numcoeff_vlc not relevant
      numcoeff_vlc = 0;	                                                                                                                                                                                                               
  
   /* subblock_x = param;
    subblock_y = param;*/
  }

  se.type  = dptype;

  se.value1 = numcoeff;
  se.value2 = numtrailingones;
  se.len    = numcoeff_vlc; /* use len to pass vlcnum */

#if TRACE
  snprintf(se.tracestring,
    TRACESTRING_SIZE, "%s # c & tr.1s(%d,%d) vlc=%d #c=%d #t1=%d",
    type, subblock_x, subblock_y, numcoeff_vlc, numcoeff, numtrailingones);
#endif

  if (!cdc)
    writeSyntaxElement_NumCoeffTrailingOnes(&se, rw_bs, block_type, bsBuff);
  else
    writeSyntaxElement_NumCoeffTrailingOnesChromaDC(p_Vid, &se, rw_bs);

  /*if (!cdc)
    writeSyntaxElement_NumCoeffTrailingOnes(&se, dataPart);
  else
    writeSyntaxElement_NumCoeffTrailingOnesChromaDC(p_Vid, &se, dataPart);*/

 // *mb_bits_coeff += se.len;
  no_bits                += se.len;

  if (!numcoeff)
    return no_bits;

  if (numcoeff)
  {
    code = 0;
    for (k = lastcoeff; k > lastcoeff - numtrailingones; k--)
    {
      level = pLevel[k]; // level
#ifdef  _DEBUG
      if (iabs(level) > 1)
      {
        printf("ERROR: level > 1\n");
        exit(-1);
      }
#endif
      code <<= 1;

      code |= (level < 0);
    }

    if (numtrailingones)
    {
      se.type  = dptype;

      se.value2 = numtrailingones;
      se.value1 = code;

#if TRACE
      snprintf(se.tracestring,
        TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)",
        type, subblock_x, subblock_y);
#endif

      writeSyntaxElement_VLC (&se, rw_bs, block_type, bsBuff);
     // writeSyntaxElement_VLC (&se, dataPart);
   //   *mb_bits_coeff += se.len;
      no_bits                += se.len;

    }

    // encode levels
    level_two_or_higher = (numcoeff > 3 && numtrailingones == 3) ? 0 : 1;

    vlcnum = (numcoeff > 10 && numtrailingones < 3) ? 1 : 0;

    for (k = lastcoeff - numtrailingones; k >= 0; k--)
    {
      level = pLevel[k]; // level

      se.value1 = level;
      se.type  = dptype;

#if TRACE
      snprintf(se.tracestring,
        TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d lev=%3d",
        type, subblock_x, subblock_y, k, vlcnum, level);
#endif

      if (level_two_or_higher)
      {
        level_two_or_higher = 0;

        if (se.value1 > 0)
          se.value1 --;
        else
          se.value1 ++;        
      }

      //    encode level

	  if (vlcnum == 0)
        writeSyntaxElement_Level_VLC1(&se, rw_bs, p_Vid->active_sps->profile_idc, block_type, bsBuff);
      else
        writeSyntaxElement_Level_VLCN(&se, vlcnum, rw_bs, p_Vid->active_sps->profile_idc, block_type, bsBuff);
      /*if (vlcnum == 0)
        writeSyntaxElement_Level_VLC1(&se, dataPart, p_Vid->active_sps->profile_idc);
      else
        writeSyntaxElement_Level_VLCN(&se, vlcnum, dataPart, p_Vid->active_sps->profile_idc);*/

      // update VLC table
      if (iabs(level) > incVlc[vlcnum])
        vlcnum++;

      if ((k == lastcoeff - numtrailingones) && iabs(level) > 3)
        vlcnum = 2;

   //   *mb_bits_coeff += se.len;
      no_bits                += se.len;
    }

    // encode total zeroes
    if (numcoeff < max_coeff_num)
    {

      se.type  = dptype;
      se.value1 = totzeros;

      vlcnum = numcoeff - 1;

      se.len = vlcnum;

#if TRACE
      snprintf(se.tracestring,
        TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d totzeros=%3d",
        type, subblock_x, subblock_y, vlcnum, totzeros);
#endif
      if (!cdc)
        writeSyntaxElement_TotalZeros(&se, rw_bs, block_type, bsBuff);
      else
        writeSyntaxElement_TotalZerosChromaDC(p_Vid, &se, rw_bs);
	  /* if (!cdc)
        writeSyntaxElement_TotalZeros(&se, dataPart);
      else
        writeSyntaxElement_TotalZerosChromaDC(p_Vid, &se, dataPart);*/

      //*mb_bits_coeff += se.len;
      no_bits                += se.len;
    }

    // encode run before each coefficient
    zerosleft = totzeros;
    numcoef = numcoeff;
    for (k = lastcoeff; k >= 0; k--)
    {
      run = pRun[k]; // run

      se.value1 = run;
      se.type   = dptype;

      // for last coeff, run is remaining totzeros
      // when zerosleft is zero, remaining coeffs have 0 run
      if ((!zerosleft) || (numcoeff <= 1 ))
        break;

      if (numcoef > 1 && zerosleft)
      {
        vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);
        se.len = vlcnum;

#if TRACE
        snprintf(se.tracestring,
          TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d run=%2d",
          type, subblock_x, subblock_y, k, vlcnum, run);
#endif

        writeSyntaxElement_Run(&se, rw_bs, block_type, bsBuff);
       // writeSyntaxElement_Run(&se, dataPart);

       // *mb_bits_coeff += se.len;
        no_bits                += se.len;

        zerosleft -= run;
        numcoef --;
      }
    }
  }

  return no_bits;
}


/*!
 ************************************************************************
 * \brief
 *    Makes code word and passes it back
 *
 * \par Input:
 *    Info   : Xn..X2 X1 X0                                             \n
 *    Length : Total number of bits in the codeword
 ************************************************************************
 */

int symbol2vlc(SyntaxElement *sym)
{
  int info_len = sym->len;

  // Convert info into a bitpattern int
  sym->bitpattern = 0;

  // vlc coding
  while(--info_len >= 0)
  {
    sym->bitpattern <<= 1;
    sym->bitpattern |= (0x01 & (sym->inf >> info_len));
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    writes UVLC code to the appropriate buffer
 ************************************************************************
 */
void  writeUVLC2buffer(SyntaxElement *se, RW_Bitstream *bs)
{
  unsigned int mask = 1 << (se->len - 1);
  byte *byte_buf  = &bs->byte_buf;
  int *bits_to_go = &bs->bits_to_go;
  int i;

  // Add the new bits to the bitstream.
  // Write out a byte if it is full
  if ( se->len < 33 )
  {
    for (i = 0; i < se->len; i++)
    {
      *byte_buf <<= 1;

      if (se->bitpattern & mask)
        *byte_buf |= 1;

      mask >>= 1;

      if ((--(*bits_to_go)) == 0)
      {
        *bits_to_go = 8;      
        bs->streamBuffer[bs->byte_pos++] = *byte_buf;
        *byte_buf = 0;      
      }
    }
  }
  else
  {
    // zeros
    for (i = 0; i < (se->len - 32); i++)
    {
      *byte_buf <<= 1;

      if ((--(*bits_to_go)) == 0)
      {
        *bits_to_go = 8;      
        bs->streamBuffer[bs->byte_pos++] = *byte_buf;
        *byte_buf = 0;      
      }
    }
    // actual info
    mask = (unsigned int) 1 << 31;
    for (i = 0; i < 32; i++)
    {
      *byte_buf <<= 1;

      if (se->bitpattern & mask)
        *byte_buf |= 1;

      mask >>= 1;

      if ((--(*bits_to_go)) == 0)
      {
        *bits_to_go = 8;      
        bs->streamBuffer[bs->byte_pos++] = *byte_buf;
        *byte_buf = 0;      
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Makes code word and passes it back
 *    A code word has the following format: 0 0 0 ... 1 Xn ...X2 X1 X0.
 *
 * \par Input:
 *    Info   : Xn..X2 X1 X0                                             \n
 *    Length : Total number of bits in the codeword
 ************************************************************************
 */
 // NOTE this function is called with sym->inf > (1<<(sym->len >> 1)).  The upper bits of inf are junk
int symbol2uvlc(SyntaxElement *sym)
{
  int suffix_len = sym->len >> 1;
  //assert (suffix_len < 32);
  suffix_len = (1 << suffix_len);
  sym->bitpattern = suffix_len | (sym->inf & (suffix_len - 1));
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    mapping for se(v) syntax elements
 * \param se
 *    value to be mapped
 * \param dummy
 *    dummy parameter
 * \param len
 *    returns mapped value length
 * \param info
 *    returns mapped value
 ************************************************************************
 */
void se_linfo(int se, int dummy, int *len,int *info)
{  
  int sign = (se <= 0) ? 1 : 0;
  int n = iabs(se) << 1;   //  n+1 is the number in the code table.  Based on this we find length and info
  int nn = (n >> 1);
  int i;
  for (i = 0; i < 33 && nn != 0; i++)
  {
    nn >>= 1;
  }
  *len  = (i << 1) + 1;
  *info = n - (1 << i) + sign;
}

/*!
 ************************************************************************
 * \brief
 *    mapping for ue(v) syntax elements
 * \param ue
 *    value to be mapped
 * \param dummy
 *    dummy parameter
 * \param info
 *    returns mapped value
 * \param len
 *    returns mapped value length
 ************************************************************************
 */
void ue_linfo(int ue, int dummy, int *len,int *info)
{
  int i, nn =(ue + 1) >> 1;

  for (i=0; i < 33 && nn != 0; i++)
  {
    nn >>= 1;
  }
  *len  = (i << 1) + 1;
  *info = ue + 1 - (1 << i);
}

void cbp_linfo_normal_intra(int cbp, int dummy, int *len,int *info)
{
  ue_linfo(RE_NCBP[1][cbp][0], dummy, len, info);
}

void cbp_linfo_normal_inter(int cbp, int dummy, int *len,int *info)
{
  ue_linfo(RE_NCBP[1][cbp][1], dummy, len, info);
}

/*!
************************************************************************
* \brief
*    generates UVLC code and passes the codeword to the buffer
* \author
*  Tian Dong
************************************************************************
*/
void writeSE_Flag(SyntaxElement *se, RW_Bitstream *rw_bs)
{
  se->len        = 1;
  se->bitpattern = (se->value1 & 1);

  writeUVLC2buffer(se, rw_bs );

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif
}


/*!
************************************************************************
* \brief
*    generates UVLC code and passes the codeword to the buffer
************************************************************************
*/
void writeCBP_VLC (Macroblock* currMB, SyntaxElement *se, RW_Bitstream *rw_bs)
{
  if (currMB->mb_type == I4MB || currMB->mb_type == SI4MB ||  currMB->mb_type == I8MB)
  {
    cbp_linfo_normal_intra (se->value1, se->value2, &(se->len), &(se->inf));
  }
  else
  {
    cbp_linfo_normal_inter (se->value1,se->value2,&(se->len),&(se->inf));
  }
  symbol2uvlc(se);

  writeUVLC2buffer(se, rw_bs);

  if(se->type != SE_HEADER)
   rw_bs->write_flag = 1;

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif
}

/*!
************************************************************************
* \brief
*    generates UVLC code and passes the codeword to the buffer
************************************************************************
*/
void writeSE_SVLC(SyntaxElement *se, RW_Bitstream *rw_bs)
{
  se_linfo (se->value1,se->value2,&(se->len),&(se->inf));
  symbol2uvlc(se);

  writeUVLC2buffer(se, rw_bs);

  if(se->type != SE_HEADER)
    rw_bs->write_flag = 1;

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif
}


/*!
************************************************************************
* \brief
*    generates UVLC code and passes the codeword to the buffer
************************************************************************
*/
void writeSE_UVLC(SyntaxElement *se, RW_Bitstream *rw_bs)
{
  ue_linfo (se->value1,se->value2,&(se->len),&(se->inf));
  symbol2uvlc(se);

  writeUVLC2buffer(se, rw_bs);

  if(se->type != SE_HEADER)
    rw_bs->write_flag = 1;

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    generates code and passes the codeword to the buffer
 ************************************************************************
 */
void writeIntraPredMode_CAVLC(SyntaxElement *se, RW_Bitstream *rw_bs)
{

  if (se->value1 == -1)
  {
    se->len = 1;
    se->inf = 1;
  }
  else
  {
    se->len = 4;
    se->inf = se->value1;
  }

  se->bitpattern = se->inf;
  writeUVLC2buffer(se, rw_bs);

  if(se->type != SE_HEADER)
    rw_bs->write_flag = 1;

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif

  return;
}


/*!
 ************************************************************************
 * \brief
 *    generates UVLC code and passes the codeword to the buffer
 * \author
 *  Tian Dong
 ************************************************************************
 */
int writeSyntaxElement2Buf_UVLC(SyntaxElement *se, RW_Bitstream *rw_bs )
{

  se->mapping(se->value1,se->value2,&(se->len),&(se->inf));

  symbol2uvlc(se);

  writeUVLC2buffer(se, rw_bs );

#if TRACE
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}



/*!
 ************************************************************************
 * \brief
 *    generates VLC code and passes the codeword to the buffer
 ************************************************************************
 */
int writeSyntaxElement_VLC(SyntaxElement *se, RW_Bitstream *rw_bs, int blocktype, RW_Bitstream *bsBuff)
{
  se->inf = se->value1;
  se->len = se->value2;
  symbol2vlc(se);

  //writeUVLC2buffer(se, rw_bs);
  ///TEST
	//writeUVLC2buffer(se, rw_bs);
  if(blocktype == LUMA_INTRA16x16DC)
	  writeUVLC2buffer(se, bsBuff);
  else
	  writeUVLC2buffer(se, rw_bs);


  /*if(se->type != SE_HEADER)
    rw_bs->write_flag = 1;*/
  if(se->type != SE_HEADER)
  {
	  if(blocktype == LUMA_INTRA16x16DC)
		  bsBuff->write_flag = 1;
	  else
		  rw_bs->write_flag = 1;
  }

#if TRACE
  if(rw_bs->bitstream->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}

/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring Blocks for Number of 
 *    Nonzero Coefficients
 *
 *    Luma Blocks
 ************************************************************************
 */
static int Repredict_nnz(Macroblock *currMB, int block_type, int i,int j)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;

  PixelPos pix;

  int pred_nnz = 0;
  int cnt      = 0;

  // left block
  get4x4Neighbour(currMB, (i << 2) - 1, (j << 2), p_Vid->mb_size[IS_LUMA], &pix);

  if ((currMB->is_intra_block == TRUE) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && (currSlice->dp_mode == PAR_DP_3))
  {
    pix.available &= currSlice->intra_block[pix.mb_addr];
    if (!pix.available)
      ++cnt;
  }

  if (pix.available)
  { 
    switch (block_type)
    {
    case LUMA:
      //pred_nnz = p_Vid->nz_coeff [pix.mb_addr ][0][pix.y][pix.x];
	  pred_nnz = rdh_nz_coeff  [pix.mb_addr ][pix.y][pix.x];
      ++cnt;
      break;
    case CB:
      pred_nnz = p_Vid->nz_coeff [pix.mb_addr ][1][pix.y][pix.x];
      ++cnt;
      break;
    case CR:
      pred_nnz = p_Vid->nz_coeff [pix.mb_addr ][2][pix.y][pix.x];
      ++cnt;
      break;
    default:
      error("writeCoeff4x4_CAVLC: Invalid block type", 600);
      break;
    }
  }

  // top block
  get4x4Neighbour(currMB, (i<<2), (j<<2) - 1, p_Vid->mb_size[IS_LUMA], &pix);

  if ((currMB->is_intra_block == TRUE) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && (currSlice->dp_mode==PAR_DP_3))
  {
    pix.available &= currSlice->intra_block[pix.mb_addr];
    if (!pix.available)
      ++cnt;
  }

  if (pix.available)
  {
    switch (block_type)
    {
    case LUMA:
      //pred_nnz += p_Vid->nz_coeff [pix.mb_addr ][0][pix.y][pix.x];
	  pred_nnz += rdh_nz_coeff  [pix.mb_addr ][pix.y][pix.x];
      ++cnt;
      break;
    case CB:
      pred_nnz += p_Vid->nz_coeff [pix.mb_addr ][1][pix.y][pix.x];
      ++cnt;
      break;
    case CR:
      pred_nnz += p_Vid->nz_coeff [pix.mb_addr ][2][pix.y][pix.x];
      ++cnt;
      break;
    default:
      error("writeCoeff4x4_CAVLC: Invalid block type", 600);
      break;
    }
  }

  if (cnt==2)
  {
    ++pred_nnz;
    pred_nnz >>= 1;
  }

  return pred_nnz;
}


/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring Blocks for Number of 
 *    Nonzero Coefficients
 *
 *    Chroma Blocks
 ************************************************************************
 */
static int Repredict_nnz_chroma(Macroblock *currMB, int i,int j)
{
  StorablePicture *dec_picture = currMB->p_Slice->dec_picture;

  if (dec_picture->chroma_format_idc != YUV444)
  {
    VideoParameters *p_Vid = currMB->p_Vid;    
    Slice *currSlice = currMB->p_Slice;
    PixelPos pix;
    int pred_nnz = 0;
    int cnt      = 0;

    //YUV420 and YUV422
    // left block
    get4x4Neighbour(currMB, ((i & 0x01)<<2) - 1, j, p_Vid->mb_size[IS_CHROMA], &pix);
	//get4x4Neighbour(currMB, ((i&0x01)<<2) - 1, j, p_Vid->mb_size[IS_CHROMA], &pix);
    if ((currMB->is_intra_block == TRUE) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && (currSlice->dp_mode==PAR_DP_3))
    {
      pix.available &= currSlice->intra_block[pix.mb_addr];
      if (!pix.available)
        ++cnt;
    }

    if (pix.available)
    {
		//pred_nnz = p_Vid->nz_coeff [pix.mb_addr ][2 * (i >> 1) + pix.x][4 + pix.y];
      pred_nnz = p_Vid->nz_coeff [pix.mb_addr ][1][pix.y][2 * (i>>1) + pix.x];
      ++cnt;
    }

    // top block
    get4x4Neighbour(currMB, ((i & 0x01)<<2), j -1, p_Vid->mb_size[IS_CHROMA],  &pix);

    if ((currMB->is_intra_block == TRUE) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && (currSlice->dp_mode==PAR_DP_3))
    {
      pix.available &= currSlice->intra_block[pix.mb_addr];
      if (!pix.available)
        ++cnt;
    }

    if (pix.available)
    {
		//pred_nnz += p_Vid->nz_coeff [pix.mb_addr ][2 * (i >> 1) + pix.x][4 + pix.y];
      pred_nnz += p_Vid->nz_coeff [pix.mb_addr ][1][pix.y][2 * (i>>1) + pix.x];
      ++cnt;
    }

    if (cnt==2)
    {
      ++pred_nnz;
      pred_nnz >>= 1;
    }
    return pred_nnz;
  }
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for NumCoeff and TrailingOnes
 ************************************************************************
 */

int writeSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *se, RW_Bitstream *rw_bs, int blocktype, RW_Bitstream *bsBuff)
{
  static const byte lentab[3][4][17] =
  {
    {   // 0702
      { 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
      { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
      { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
      { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16},
    },
    {
      { 2, 6, 6, 7, 8, 8, 9,11,11,12,12,12,13,13,13,14,14},
      { 0, 2, 5, 6, 6, 7, 8, 9,11,11,12,12,13,13,14,14,14},
      { 0, 0, 3, 6, 6, 7, 8, 9,11,11,12,12,13,13,13,14,14},
      { 0, 0, 0, 4, 4, 5, 6, 6, 7, 9,11,11,12,13,13,13,14},
    },
    {
      { 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9,10,10,10,10},
      { 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10},
      { 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,10},
      { 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9,10,10,10},
    },

  };

  static const byte codtab[3][4][17] =
  {
    {
      { 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7,4},
      { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10,6},
      { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9,5},
      { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12,8},
    },
    {
      { 3,11, 7, 7, 7, 4, 7,15,11,15,11, 8,15,11, 7, 9,7},
      { 0, 2, 7,10, 6, 6, 6, 6,14,10,14,10,14,10,11, 8,6},
      { 0, 0, 3, 9, 5, 5, 5, 5,13, 9,13, 9,13, 9, 6,10,5},
      { 0, 0, 0, 5, 4, 6, 8, 4, 4, 4,12, 8,12,12, 8, 1,4},
    },
    {
      {15,15,11, 8,15,11, 9, 8,15,11,15,11, 8,13, 9, 5,1},
      { 0,14,15,12,10, 8,14,10,14,14,10,14,10, 7,12, 8,4},
      { 0, 0,13,14,11, 9,13, 9,13,10,13, 9,13, 9,11, 7,3},
      { 0, 0, 0,12,11,10, 9, 8,13,12,12,12, 8,12,10, 6,2},
    },
  };
  int vlcnum = se->len;

  // se->value1 : numcoeff
  // se->value2 : numtrailingones

  if (vlcnum == 3)
  {
    se->len = 6;  // 4 + 2 bit FLC
    if (se->value1 > 0)
    {
      se->inf = ((se->value1-1) << 2) | se->value2;
    }
    else
    {
      se->inf = 3;
    }
  }
  else
  {
    se->len = lentab[vlcnum][se->value2][se->value1];
    se->inf = codtab[vlcnum][se->value2][se->value1];
  }

  if (se->len == 0)
  {
    printf("ERROR: (numcoeff,trailingones) not valid: vlc=%d (%d, %d)\n",
      vlcnum, se->value1, se->value2);
    exit(-1);
  }

  symbol2vlc(se);

///TEST
	//writeUVLC2buffer(se, rw_bs);
  if(blocktype == LUMA_INTRA16x16DC)
	  writeUVLC2buffer(se, bsBuff);
  else
	  writeUVLC2buffer(se, rw_bs);

  

  if(se->type != SE_HEADER)
  {
	  if(blocktype == LUMA_INTRA16x16DC)
		  bsBuff->write_flag = 1;
	  else
		  rw_bs->write_flag = 1;
  }

#if TRACE
  if(rw_bs->bitstream->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for NumCoeff and TrailingOnes for Chroma DC
 ************************************************************************
 */
int writeSyntaxElement_NumCoeffTrailingOnesChromaDC(VideoParameters *p_Vid, SyntaxElement *se, RW_Bitstream *rw_bs)
{
  static const byte lentab[3][4][17] =
  {
    //YUV420
   {{ 2, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 1, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 3, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV422
   {{ 1, 7, 7, 9, 9,10,11,12,13, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 2, 7, 7, 9,10,11,12,12, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 3, 7, 7, 9,10,11,12, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 5, 6, 7, 7,10,11, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV444
   {{ 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
    { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
    { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
    { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16}}
  };

  static const byte codtab[3][4][17] =
  {
    //YUV420
   {{ 1, 7, 4, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 1, 6, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV422
   {{ 1,15,14, 7, 6, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 1,13,12, 5, 6, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 1,11,10, 4, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 1, 1, 9, 8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0}},
    //YUV444
   {{ 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7, 4},
    { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10, 6},
    { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9, 5},
    { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12, 8}}

  };  
  int yuv = p_Vid->yuv_format - 1;

  // se->value1 : numcoeff
  // se->value2 : numtrailingones
  se->len = lentab[yuv][se->value2][se->value1];
  se->inf = codtab[yuv][se->value2][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (numcoeff,trailingones) not valid: (%d, %d)\n",
      se->value1, se->value2);
    exit(-1);
  }

  symbol2vlc(se);

  writeUVLC2buffer(se, rw_bs);

  if(se->type != SE_HEADER)
    rw_bs->write_flag = 1;

#if TRACE
  if(rw_bs->bitstream->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for TotalZeros
 ************************************************************************
 */
int writeSyntaxElement_TotalZeros(SyntaxElement *se, RW_Bitstream *rw_bs, int blocktype, RW_Bitstream *bsBuff)
{
  static const byte lentab[TOTRUN_NUM][16] =
  {
    { 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},
    { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},
    { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},
    { 5,3,4,4,3,3,3,4,3,4,5,5,5},
    { 4,4,4,3,3,3,3,3,4,5,4,5},
    { 6,5,3,3,3,3,3,3,4,3,6},
    { 6,5,3,3,3,2,3,4,3,6},
    { 6,4,5,3,2,2,3,3,6},
    { 6,6,4,2,2,3,2,5},
    { 5,5,3,2,2,2,4},
    { 4,4,3,3,1,3},
    { 4,4,2,1,3},
    { 3,3,1,2},
    { 2,2,1},
    { 1,1},
  };

  static const byte codtab[TOTRUN_NUM][16] =
  {
    {1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
    {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
    {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
    {3,7,5,4,6,5,4,3,3,2,2,1,0},
    {5,4,3,7,6,5,4,3,2,1,1,0},
    {1,1,7,6,5,4,3,2,1,1,0},
    {1,1,5,4,3,3,2,1,1,0},
    {1,1,1,3,3,2,2,1,0},
    {1,0,1,3,2,1,1,1,},
    {1,0,1,3,2,1,1,},
    {0,1,1,2,1,3},
    {0,1,1,1,1},
    {0,1,1,1},
    {0,1,1},
    {0,1},
  };
  int vlcnum = se->len;

  // se->value1 : TotalZeros
  se->len = lentab[vlcnum][se->value1];
  se->inf = codtab[vlcnum][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (TotalZeros) not valid: (%d)\n",se->value1);
    exit(-1);
  }

  symbol2vlc(se);

  ///TEST
	//writeUVLC2buffer(se, rw_bs);
  if(blocktype == LUMA_INTRA16x16DC)
	  writeUVLC2buffer(se, bsBuff);
  else
	  writeUVLC2buffer(se, rw_bs);


  if(se->type != SE_HEADER)
  {
	  if(blocktype == LUMA_INTRA16x16DC)
		  bsBuff->write_flag = 1;
	  else
		  rw_bs->write_flag = 1;
  }

#if TRACE
  if(rw_bs->bitstream->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for TotalZeros for Chroma DC
 ************************************************************************
 */
int writeSyntaxElement_TotalZerosChromaDC(VideoParameters *p_Vid, SyntaxElement *se, RW_Bitstream *rw_bs)
{
  static const byte lentab[3][TOTRUN_NUM][16] =
  {
    //YUV420
   {{ 1,2,3,3},
    { 1,2,2},
    { 1,1}},
    //YUV422
   {{ 1,3,3,4,4,4,5,5},
    { 3,2,3,3,3,3,3},
    { 3,3,2,2,3,3},
    { 3,2,2,2,3},
    { 2,2,2,2},
    { 2,2,1},
    { 1,1}},
    //YUV444
   {{ 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},
    { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},
    { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},
    { 5,3,4,4,3,3,3,4,3,4,5,5,5},
    { 4,4,4,3,3,3,3,3,4,5,4,5},
    { 6,5,3,3,3,3,3,3,4,3,6},
    { 6,5,3,3,3,2,3,4,3,6},
    { 6,4,5,3,2,2,3,3,6},
    { 6,6,4,2,2,3,2,5},
    { 5,5,3,2,2,2,4},
    { 4,4,3,3,1,3},
    { 4,4,2,1,3},
    { 3,3,1,2},
    { 2,2,1},
    { 1,1}}
  };

  static const byte codtab[3][TOTRUN_NUM][16] =
  {
    //YUV420
   {{ 1,1,1,0},
    { 1,1,0},
    { 1,0}},
    //YUV422
   {{ 1,2,3,2,3,1,1,0},
    { 0,1,1,4,5,6,7},
    { 0,1,1,2,6,7},
    { 6,0,1,2,7},
    { 0,1,2,3},
    { 0,1,1},
    { 0,1}},
    //YUV444
   {{1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
    {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
    {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
    {3,7,5,4,6,5,4,3,3,2,2,1,0},
    {5,4,3,7,6,5,4,3,2,1,1,0},
    {1,1,7,6,5,4,3,2,1,1,0},
    {1,1,5,4,3,3,2,1,1,0},
    {1,1,1,3,3,2,2,1,0},
    {1,0,1,3,2,1,1,1,},
    {1,0,1,3,2,1,1,},
    {0,1,1,2,1,3},
    {0,1,1,1,1},
    {0,1,1,1},
    {0,1,1},
    {0,1}}
  };
  int vlcnum = se->len;
  int yuv = p_Vid->yuv_format - 1;  

  // se->value1 : TotalZeros
  se->len = lentab[yuv][vlcnum][se->value1];
  se->inf = codtab[yuv][vlcnum][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (TotalZeros) not valid: (%d)\n",se->value1);
    exit(-1);
  }

  symbol2vlc(se);

  writeUVLC2buffer(se, rw_bs);

  if(se->type != SE_HEADER)
    rw_bs->write_flag = 1;

#if TRACE
  if(rw_bs->bitstream->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for Run Before Next Coefficient, VLC0
 ************************************************************************
 */
int writeSyntaxElement_Run(SyntaxElement *se, RW_Bitstream *rw_bs, int blocktype, RW_Bitstream *bsBuff)
{
  static const byte lentab[TOTRUN_NUM][16] =
  {
    {1,1},
    {1,2,2},
    {2,2,2,2},
    {2,2,2,3,3},
    {2,2,3,3,3,3},
    {2,3,3,3,3,3,3},
    {3,3,3,3,3,3,3,4,5,6,7,8,9,10,11},
  };

  static const byte codtab[TOTRUN_NUM][16] =
  {
    {1,0},
    {1,1,0},
    {3,2,1,0},
    {3,2,1,1,0},
    {3,2,3,2,1,0},
    {3,0,1,3,2,5,4},
    {7,6,5,4,3,2,1,1,1,1,1,1,1,1,1},
  };
  int vlcnum = se->len;

  // se->value1 : run
  se->len = lentab[vlcnum][se->value1];
  se->inf = codtab[vlcnum][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (run) not valid: (%d)\n",se->value1);
    exit(-1);
  }

  symbol2vlc(se);

  ///TEST
	//writeUVLC2buffer(se, rw_bs);
  if(blocktype == LUMA_INTRA16x16DC)
	  writeUVLC2buffer(se, bsBuff);
  else
	  writeUVLC2buffer(se, rw_bs);


  if(se->type != SE_HEADER)
  {
	  if(blocktype == LUMA_INTRA16x16DC)
		  bsBuff->write_flag = 1;
	  else
		  rw_bs->write_flag = 1;
  }

#if TRACE
  if(rw_bs->bitstream->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for Coeff Level (VLC1)
 ************************************************************************
 */
int writeSyntaxElement_Level_VLC1(SyntaxElement *se, RW_Bitstream *rw_bs, int profile_idc, int blocktype, RW_Bitstream *bsBuff)
{
  int level  = se->value1;
  int sign   = (level < 0 ? 1 : 0);
  int levabs = iabs(level);

  if (levabs < 8)
  {
    se->len = levabs * 2 + sign - 1;
    se->inf = 1;
  }
  else if (levabs < 16) 
  {
    // escape code1
    se->len = 19;
    se->inf = 16 | ((levabs << 1) - 16) | sign;
  }
  else
  {
    int iMask = 4096, numPrefix = 0;
    int levabsm16 = levabs + 2032;

    // escape code2
    if ((levabsm16) >= 4096)
    {
      numPrefix++;
      while ((levabsm16) >= (4096 << numPrefix))
      {
        numPrefix++;
      }
    }
   
    iMask <<= numPrefix;
    se->inf = iMask | ((levabsm16 << 1) - iMask) | sign;

    /* Assert to make sure that the code fits in the VLC */
    /* make sure that we are in High Profile to represent level_prefix > 15 */
    if (numPrefix > 0 && !is_FREXT_profile( profile_idc ))
    {
      //error( "level_prefix must be <= 15 except in High Profile\n",  1000 );
      se->len = 0x0000FFFF; // This can be some other big number
      return (se->len);
    }
    
    se->len = 28 + (numPrefix << 1);
  }

  symbol2vlc(se);

  //writeUVLC2buffer(se, rw_bs);
  ///TEST
	//writeUVLC2buffer(se, rw_bs);
  if(blocktype == LUMA_INTRA16x16DC)
	  writeUVLC2buffer(se, bsBuff);
  else
	  writeUVLC2buffer(se, rw_bs);


  if(se->type != SE_HEADER)
  {
	  if(blocktype == LUMA_INTRA16x16DC)
		  bsBuff->write_flag = 1;
	  else
		  rw_bs->write_flag = 1;
  }

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for Coeff Level
 ************************************************************************
 */
int writeSyntaxElement_Level_VLCN(SyntaxElement *se, int vlc, RW_Bitstream *rw_bs, int profile_idc, int blocktype, RW_Bitstream *bsBuff)
{  
  int level  = se->value1;
  int sign   = (level < 0 ? 1 : 0);
  int levabs = iabs(level) - 1;  

  int shift = vlc - 1;        
  int escape = (15 << shift);

  if (levabs < escape)
  {
    int sufmask   = ~((0xffffffff) << shift);
    int suffix    = (levabs) & sufmask;

    se->len = ((levabs) >> shift) + 1 + vlc;
    se->inf = (2 << shift) | (suffix << 1) | sign;
  }
  else
  {
    int iMask = 4096;
    int levabsesc = levabs - escape + 2048;
    int numPrefix = 0;

    if ((levabsesc) >= 4096)
    {
      numPrefix++;
      while ((levabsesc) >= (4096 << numPrefix))
      {
        numPrefix++;
      }
    }

    iMask <<= numPrefix;
    se->inf = iMask | ((levabsesc << 1) - iMask) | sign;

    /* Assert to make sure that the code fits in the VLC */
    /* make sure that we are in High Profile to represent level_prefix > 15 */
    if (numPrefix > 0 &&  !is_FREXT_profile( profile_idc ))
    {
      //error( "level_prefix must be <= 15 except in High Profile\n",  1000 );
      se->len = 0x0000FFFF; // This can be some other big number
      return (se->len);
    }
    se->len = 28 + (numPrefix << 1);
  }

  symbol2vlc(se);

  ///TEST
	//writeUVLC2buffer(se, rw_bs);
  if(blocktype == LUMA_INTRA16x16DC)
	  writeUVLC2buffer(se, bsBuff);
  else
	  writeUVLC2buffer(se, rw_bs);


  if(se->type != SE_HEADER)
  {
	  if(blocktype == LUMA_INTRA16x16DC)
		  bsBuff->write_flag = 1;
	  else
		  rw_bs->write_flag = 1;
  }

#if TRACE
  if(rw_bs->trace_enabled)
    trace2out (se);
#endif

  return (se->len);
}


 /*!
 ************************************************************************
 * \brief
 *    Converts String Of Data Bits (SODB) to Raw Byte Sequence
 *    Packet (RBSP)
 * \param currStream
 *        Bitstream which contains data bits.
 * \return None
 * \note currStream is byte-aligned at the end of this function
 *
 ************************************************************************
*/

void SODBtoRBSP(RW_Bitstream *rw_bs)
{
	//if(rw_bs->bits_to_go !=8)
	{
		rw_bs->byte_buf <<= 1;
		rw_bs->byte_buf |= 1;
		rw_bs->bits_to_go--;
		rw_bs->byte_buf <<= rw_bs->bits_to_go;
		rw_bs->streamBuffer[rw_bs->byte_pos++] = rw_bs->byte_buf;
		rw_bs->bits_to_go = 8;
		rw_bs->byte_buf = 0;
	}
}



int scie_rdh_4x4(int****my_cof, int block_x4, int block_y4, int mb_num)  
{
	int numcoeff=0;
	int x,y,i,j;

	int my_coef[2];
	char sec[4];
	int k=0;
	int numk;

	i = mb_num;
	x = block_x4;
	y = block_y4;

	if(my_cof[i][PLANE_Y][x+1][y+3]==0 && my_cof[i][PLANE_Y][x+2][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+1]==0 && my_cof[i][PLANE_Y][x+2][y+3]==0 && my_cof[i][PLANE_Y][x+3][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+3]==0)
	{
#if EMB_CTRL
		int k = 0;
		//				int numk;
		//*/
		int RDHok=0;

		switch(DHrate)
		{
		case 6:
			RDHok=1;
			break;
		case 7:
			if(my_cof[i][PLANE_Y][x][y+2]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		case 8:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		case 9:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		case 10:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y]==0 && my_cof[i][PLANE_Y][x+1][y]==0)
				RDHok=1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
			else
				RDHok=0;
			break;
		case 11:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y]==0 && my_cof[i][PLANE_Y][x+1][y]==0 && my_cof[i][PLANE_Y][x][y+1]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		}
		if(DHpairing == 4 && RDHok==1)
#endif				
		{
			//디코딩
			//for(j=3;j>-1;j--)
			//인코딩
			int start, end;

			if(rdh_encode == 1) { start = 0; end= 4;}
			else if (rdh_encode == 2) { start = 3; end= -1;}

			j= start;
			while(j!=end)
			{
				int dh;
				int x1,y1,x2,y2;

				x1 = x + j;
				x2 = x + ((j+1)%4);
				y1 = y + (3-j);
				y2 = y + (3-((j+1)%4));


				dh = (int)ceil(((double)(my_cof[i][PLANE_Y][x1][y1] - (double)my_cof[i][PLANE_Y][x2][y2])/(double)2));
				MB_emb_flag = 1;
				//#if RDH_ENCOD
				if(rdh_encode == 1) {
					if(dh == 0)
					{
						char secbit;

						if(randnum[rand_count] == '1')
							secbit = 1;
						else
							secbit = 0;

						my_cof[i][PLANE_Y][x1][y1] = my_cof[i][PLANE_Y][x1][y1] + secbit;
						my_cof[i][PLANE_Y][x2][y2] = my_cof[i][PLANE_Y][x2][y2] - secbit;



						rand_count++;
						rdh_count++;
					}
					else if(dh > 0)
					{
						my_cof[i][PLANE_Y][x1][y1] = my_cof[i][PLANE_Y][x1][y1] + 1;
						my_cof[i][PLANE_Y][x2][y2] = my_cof[i][PLANE_Y][x2][y2] - 1;
					}
				}
				//#else
				else if(rdh_encode == 2) {
					if(dh == 0)
					{
						sec[k] = 0;
						//fprintf(fp_sec,"%d",sec);
						num_sec++;

						MB_emb_flag = 1;
						my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
						my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];
						k++;
					}
					else if(dh == 1)
					{
						sec[k] = 1;
						//fprintf(fp_sec,"%d",sec);
						num_sec++;
						MB_emb_flag = 1;
						my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
						my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];

						my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
						my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
						k++;
					}
					else if(dh > 1)
					{
						my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - 1;
						my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + 1;

						my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
						my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
					}
					else if(dh < 0)
					{
						my_coef[0] = my_cof[i][PLANE_Y][x1][y1];
						my_coef[1] = my_cof[i][PLANE_Y][x2][y2];
					}
				}
				//#endif
				if(rdh_encode == 1) j++;
				else if (rdh_encode == 2) j--;
			}
			//#if RDH_ENCOD
			//				for(j=0;j<4;j++)
			//#else
			//				for(j=3;j>-1;j--)
			//#endif
			//				{
			//					int dh;
			//					int x1,y1,x2,y2;
			//
			//					x1 = x + j;
			//					x2 = x + ((j+1)%4);
			//					y1 = y + (3-j);
			//					y2 = y + (3-((j+1)%4));
			//
			//
			//					dh = (int)ceil(((double)(my_cof[i][PLANE_Y][x1][y1] - (double)my_cof[i][PLANE_Y][x2][y2])/(double)2));
			//
			//#if RDH_ENCOD
			//					if(dh == 0)
			//					{
			//						char secbit;
			//
			//						if(randnum[rand_count] == '1')
			//							secbit = 1;
			//						else
			//							secbit = 0;
			//
			//						my_cof[i][PLANE_Y][x1][y1] = my_cof[i][PLANE_Y][x1][y1] + secbit;
			//						my_cof[i][PLANE_Y][x2][y2] = my_cof[i][PLANE_Y][x2][y2] - secbit;
			//
			//						rand_count++;
			//						rdh_count++;
			//					}
			//					else if(dh > 0)
			//					{
			//						my_cof[i][PLANE_Y][x1][y1] = my_cof[i][PLANE_Y][x1][y1] + 1;
			//						my_cof[i][PLANE_Y][x2][y2] = my_cof[i][PLANE_Y][x2][y2] - 1;
			//					}
			//
			//#else
			//					if(dh == 0)
			//						{
			//							sec[k] = 0;
			//							//fprintf(fp_sec,"%d",sec);
			//							num_sec++;
			//
			//							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
			//							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];
			//							k++;
			//						}
			//						else if(dh == 1)
			//						{
			//							sec[k] = 1;
			//							//fprintf(fp_sec,"%d",sec);
			//							num_sec++;
			//
			//							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
			//							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];
			//
			//							my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
			//							my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
			//							k++;
			//						}
			//						else if(dh > 1)
			//						{
			//							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - 1;
			//							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + 1;
			//
			//							my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
			//							my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
			//						}
			//						else if(dh < 0)
			//						{
			//							my_coef[0] = my_cof[i][PLANE_Y][x1][y1];
			//							my_coef[1] = my_cof[i][PLANE_Y][x2][y2];
			//						}
			//#endif
			//
			//				}



#if TEST1
			else if(DHpairing == 3)
			{
				if(RDHok==1)
				{
					for(j=2;j>-1;j--)
					{
						int dh;
						int x1,y1,x2,y2;


						x1 = x+j;
						x2 = x+((j+1)%3);
						y1 = y+(3-j);
						y2 = y+(3-((j+1)%3));
						//
#if RDH_DE_Doub_mod
						dh = (int)ceil(((double)(my_cof[i][PLANE_Y][x1][y1] - (double)my_cof[i][PLANE_Y][x2][y2])/(double)2));
#else
						dh = (int)floor(((double)(my_cof[i][PLANE_Y][x1][y1] - (double)my_cof[i][PLANE_Y][x2][y2])/(double)2));
#endif
						if(dh == 0)
						{
							sec[k] = 0;
							////fprintf(fp_sec,"%d",sec);
							//num_sec++;

							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];
							////k++;
						}
						else if(dh == 1)
						{
							sec[k] = 1;
							////fprintf(fp_sec,"%d",sec);
							//num_sec++;

							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];

							my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
							my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
							////k++;
						}
						else if(dh > 1)
						{
							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - 1;
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + 1;

							my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
							my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
						}
						else if(dh < 0)
						{
							my_coef[0] = my_cof[i][PLANE_Y][x1][y1];
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2];
						}


						// inverse quant for 4x4 transform only
						currSlice->cof[0][x1][y1]= rshift_rnd_sf((my_coef[0] * InvLevelScale4x4[x1-x][y1-y])<<qp_per, 4);
						currSlice->cof[0][x2][y2]= rshift_rnd_sf((my_coef[1] * InvLevelScale4x4[x2-x][y2-y])<<qp_per, 4);


					}
				}
			}
#endif

			//#if !RDH_ENCOD
			if(rdh_encode == 2) {
				for(numk=k-1;numk>=0;numk--)
					fprintf(fp_sec,"%d",sec[numk]);
			}
			//#endif
		}
	}		
	return numcoeff;
}



int sci_rdh_4x4(int****my_cof, int block_x4, int block_y4, int mb_num)  
{
	int numcoeff=0;
	int x,y,i,j;

	int my_coef[2];
	char sec[4];
	int k=0;
	int numk;

	i = mb_num;
	x = block_x4;
	y = block_y4;

	if(my_cof[i][PLANE_Y][x+1][y+3]==0 && my_cof[i][PLANE_Y][x+2][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+1]==0 && my_cof[i][PLANE_Y][x+2][y+3]==0 && my_cof[i][PLANE_Y][x+3][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+3]==0)
	{
#if EMB_CTRL
		int k = 0;
		//				int numk;
		//*/
		int RDHok=0;

		switch(DHrate)
		{
		case 6:
			RDHok=1;
			break;
		case 7:
			if(my_cof[i][PLANE_Y][x][y+2]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		case 8:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		case 9:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		case 10:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y]==0 && my_cof[i][PLANE_Y][x+1][y]==0)
				RDHok=1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
			else
				RDHok=0;
			break;
		case 11:
			if(my_cof[i][PLANE_Y][x][y+2]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y]==0 && my_cof[i][PLANE_Y][x+1][y]==0 && my_cof[i][PLANE_Y][x][y+1]==0)
				RDHok=1;
			else
				RDHok=0;
			break;
		}
		if(DHpairing == 4 && RDHok==1)
#endif				

#if 0//0603
		{
			int r7 = my_cof[i][PLANE_Y][x][y+3];
			int r8 = my_cof[i][PLANE_Y][x+1][y+2];
			int r9 = my_cof[i][PLANE_Y][x+2][y+1];
			int r10 = my_cof[i][PLANE_Y][x+3][y];

			totalRDok++;
			//#if RDH_ENCOD
			if( rdh_encode ==1) {
				if(r7==0 && r8==0 && r9==0 && r10==0 )
				{
					RDok[0]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') // 0000비트에 101을 숨길 때 
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1') // 0000비트에 101을 숨길 때 
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=3; rdh_count +=3;
						return numcoeff;
					}
					//if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') // 0000비트에 101을 숨길 때 
					//{
					//	my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
					//	rand_count +=4; rdh_count +=4;
					//	return numcoeff;
					//}
					//if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					//{
					//	my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
					//	rand_count +=4; rdh_count +=4;
					//	return numcoeff;
					//}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') // 0000비트에 101을 숨길 때 
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}

				}
				if(r7==-1 && r8==0 && r9==0 && r10==0 )
				{
					RDok[1]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}

				}
				if(r7==0 && r8==0 && r9==0 && r10==1 )
				{
					RDok[2]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = -1;	my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==0 && r9==1 && r10==0 )
				{
					RDok[3]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==1 && r9==0 && r10==0 )
				{
					RDok[4]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0;	my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				//if(r7==-1 && r8==-1 && r9==0 && r10==0 )
				//{
				//	//RDok[3]++;
				//	if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
				//	{
				//		my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
				//		rand_count +=4; rdh_count +=4;
				//		return numcoeff;
				//	}
				//}
				//if(r7==0 && r8==1 && r9==0 && r10==1 )
				//{
				//	//RDok[3]++;
				//	if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
				//	{
				//		my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
				//		rand_count +=4; rdh_count +=4;
				//		return numcoeff;
				//	}
				//}
				if(r7==-1 && r8==0 && r9==0 && r10==1 )
				{
					//RDok[3]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==0 && r9==1 && r10==1 )
				{
					//RDok[3]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
			}
			//#else
			else if ( rdh_encode == 2) {
				if(r7==0 && r8==1 && r9==0 && r10==-1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==-1 && r9==0 && r10==0 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==-1 && r9==1 && r10==-1 )
				{
					//RDok[0]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==0 && r9==1 && r10==0 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==1 && r9==0 && r10==0 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;

					sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==1 && r9==-1 && r10==1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==0 && r9==1 && r10==-1 )
				{
					//RDok[2]++;
					my_cof[i][PLANE_Y][x][y+3] = 0;	my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==-1 && r9==0 && r10==1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==0 && r9==-1 && r10==1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==0 && r9==-1 && r10==0 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				//////////////////////
				if((r7==0 && r8==0 && r9==0 && r10==-1) || (r7==0 && r8==0 && r9==-1 && r10==0) || (r7==0 && r8==-1 && r9==0 && r10==0) || (r7==-1 && r8==0 && r9==0 && r10==0))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==0 && r8==0 && r9==0 && r10==-1 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 1;  num_sec+=3;}

					if(r7==0 && r8==0 && r9==-1 && r10==0 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0; num_sec+=4;}

					if(r7==0 && r8==-1 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; num_sec+=4;}

					if(r7==-1 && r8==0 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; num_sec+=4;}

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if((r7==0 && r8==0 && r9==0 && r10==1 ) || (r7==0 && r8==0 && r9==1 && r10==0 ) || (r7==0 && r8==1 && r9==0 && r10==0 ) || (r7==1 && r8==0 && r9==0 && r10==0 ))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==0 && r8==0 && r9==0 && r10==1 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					if(r7==0 && r8==0 && r9==1 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; num_sec+=4;}

					if(r7==0 && r8==1 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					if(r7==1 && r8==0 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; num_sec+=4;}

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if((r7==1 && r8==0 && r9==0 && r10==1 ) || (r7==-1 && r8==0 && r9==0 && r10==-1 ))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==1 && r8==0 && r9==0 && r10==1 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; num_sec+=4;}

					if(r7==-1 && r8==0 && r9==0 && r10==-1 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
			}
			//#endif
#endif	
#if 0//160913
		{
			int r7 = my_cof[i][PLANE_Y][x][y+3];
			int r8 = my_cof[i][PLANE_Y][x+1][y+2];
			int r9 = my_cof[i][PLANE_Y][x+2][y+1];
			int r10 = my_cof[i][PLANE_Y][x+3][y];

			totalRDok++;
			//#if RDH_ENCOD
			if(rdh_encode == 1) {
				if(r7==0 && r8==0 && r9==0 && r10==0 )
				{
					RDok[0]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') // 0000비트에 101을 숨길 때 
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' ) //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=3; rdh_count +=3;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' & randnum[rand_count+3] == '0') // 0000비트에 101을 숨길 때 
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					//if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') // 0000비트에 101을 숨길 때 
					//{
					//	my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
					//	rand_count +=4; rdh_count +=4;
					//	return numcoeff;
					//}
					//if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0000 비트에  0101을 숨기면 0010으로
					//{
					//	my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
					//	rand_count +=4; rdh_count +=4;
					//	return numcoeff;
					//}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') // 0000비트에 101을 숨길 때 
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1') //0000 비트에  0101을 숨기면 0010으로
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=3; rdh_count +=3;
						return numcoeff;
					}

				}
				if(r7==-1 && r8==0 && r9==0 && r10==0 )
				{
					RDok[1]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1') //0100 비트에  1000을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}

				}
				if(r7==0 && r8==0 && r9==0 && r10==1 )
				{
					RDok[2]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0;	my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==0 && r9==1 && r10==0 )
				{
					RDok[3]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==1 && r9==0 && r10==0 )
				{
					RDok[4]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 1;	my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				//if(r7==-1 && r8==-1 && r9==0 && r10==0 )
				//{
				//	//RDok[3]++;
				//	if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
				//	{
				//		my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
				//		rand_count +=4; rdh_count +=4;
				//		return numcoeff;
				//	}
				//}
				//if(r7==0 && r8==1 && r9==0 && r10==1 )
				//{
				//	//RDok[3]++;
				//	if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
				//	{
				//		my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
				//		rand_count +=4; rdh_count +=4;
				//		return numcoeff;
				//	}
				//}
				if(r7==-1 && r8==0 && r9==-1 && r10==0 )
				{
					//RDok[3]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==1 && r9==1 && r10==0 )
				{
					//RDok[3]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0') //0001 비트에  0010을 숨기면
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
			}
			//#else
			if(rdh_encode == 2) {
				if(r7==-1 && r8==1 && r9==-1 && r10==1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==-1 && r9==0 && r10==1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==0 && r9==-1 && r10==0 )
				{
					//RDok[0]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==0 && r9==1 && r10==-1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==1 && r9==0 && r10==-1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;

					sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==0 && r9==1 && r10==0 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==0 && r9==0 && r10==1 )
				{
					//RDok[2]++;
					my_cof[i][PLANE_Y][x][y+3] = 0;	my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==-1 && r9==1 && r10==-1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==0 && r9==0 && r10==-1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==0 && r9==-1 && r10==1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				//////////////////////
				if((r7==0 && r8==0 && r9==0 && r10==-1) || (r7==0 && r8==0 && r9==-1 && r10==0) || (r7==0 && r8==-1 && r9==0 && r10==0) || (r7==-1 && r8==0 && r9==0 && r10==0))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==0 && r8==0 && r9==0 && r10==-1 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; num_sec+=4;}

					if(r7==0 && r8==0 && r9==-1 && r10==0 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					if(r7==0 && r8==-1 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; num_sec+=4;}

					if(r7==-1 && r8==0 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; num_sec+=4;}

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if((r7==0 && r8==0 && r9==0 && r10==1 ) || (r7==0 && r8==0 && r9==1 && r10==0 ) || (r7==0 && r8==1 && r9==0 && r10==0 ) || (r7==1 && r8==0 && r9==0 && r10==0 ))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==0 && r8==0 && r9==0 && r10==1 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=3;}

					if(r7==0 && r8==0 && r9==1 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; num_sec+=4;}

					if(r7==0 && r8==1 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; num_sec+=4;}

					if(r7==1 && r8==0 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if((r7==1 && r8==1 && r9==0 && r10==0 ) || (r7==-1 && r8==-1 && r9==0 && r10==0 ))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==1 && r8==1 && r9==0 && r10==0 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					if(r7==-1 && r8==-1 && r9==0 && r10==0 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; num_sec+=3;}

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
			}
			//#endif
#endif

#if 1
		{
			int r7 = my_cof[i][PLANE_Y][x][y+3];
			int r8 = my_cof[i][PLANE_Y][x+1][y+2];
			int r9 = my_cof[i][PLANE_Y][x+2][y+1];
			int r10 = my_cof[i][PLANE_Y][x+3][y];

			totalRDok++;

			//#if RDH_ENCOD
			if(rdh_encode == 1)
			{
				if(r7==0 && r8==0 && r9==0 && r10==0 )
				{
					RDok[0]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=3; rdh_count +=3;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=3; rdh_count +=3;
						return numcoeff;
					}
#if 1//PS
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
#endif			

#if 1//FS
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
#endif

				}
				if(r7==-1 && r8==0 && r9==0 && r10==0 )
				{
					RDok[1]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = -1; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}

				}
				if(r7==0 && r8==0 && r9==0 && r10==1 )
				{
					RDok[2]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0;	my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==0 && r9==1 && r10==0 )
				{
					RDok[3]++;
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
				if(r7==0 && r8==1 && r9==0 && r10==0 )
				{
					RDok[4]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = 1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0;	my_cof[i][PLANE_Y][x+3][y] = -1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
#if 0//PS				
				if(r7==0 && r8==1 && r9==1 && r10==0 )
				{
					RDok[3]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' &&  randnum[rand_count+2] == '0' && randnum[rand_count+3] == '0')
					{
						my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 1;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
#endif
#if 0//FS
				if(r7==-1 && r8==0 && r9==-1 && r10==0 )
				{
					RDok[3]++;
					if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' &&  randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1')
					{
						my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;
						rand_count +=4; rdh_count +=4;
						return numcoeff;
					}
				}
#endif
			}			
			//#else
			else if(rdh_encode == 2) 
			{
				if(r7==0 && r8==-1 && r9==0 && r10==1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==0 && r9==1 && r10==-1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==0 && r9==0 && r10==1 )
				{
					//RDok[0]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==-1 && r9==1 && r10==-1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==1 && r9==-1 && r10==0 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 1;

					sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==0 && r8==1 && r9==0 && r10==-1 )
				{
					//RDok[1]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==-1 && r8==1 && r9==-1 && r10==1 )
				{
					//RDok[2]++;
					my_cof[i][PLANE_Y][x][y+3] = 0;	my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if(r7==1 && r8==0 && r9==0 && r10==-1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
#if 0//PS
				if(r7==0 && r8==0 && r9==-1 && r10==1 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 1; my_cof[i][PLANE_Y][x+2][y+1] = 1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
#endif
#if 0//FS
				if(r7==-1 && r8==0 && r9==1 && r10==0 )
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = -1; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = -1; my_cof[i][PLANE_Y][x+3][y] = 0;

					sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; sec[k++] = 1;
					num_sec+=4;
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
#endif

				//////////////////////
				if((r7==0 && r8==0 && r9==0 && r10==-1) || (r7==0 && r8==0 && r9==-1 && r10==0) || (r7==0 && r8==-1 && r9==0 && r10==0) || (r7==-1 && r8==0 && r9==0 && r10==0) || (r7==-1 && r8==-1 && r9==0 && r10==0)) //FS
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==0 && r8==0 && r9==0 && r10==-1 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					if(r7==0 && r8==0 && r9==-1 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; num_sec+=4;}

					if(r7==0 && r8==-1 && r9==0 && r10==0 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}

					if(r7==-1 && r8==0 && r9==0 && r10==0 )
					{sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=3;}
#if 1//FS
					if(r7==-1 && r8==-1 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; num_sec+=4;}
#endif

					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
				if((r7==0 && r8==0 && r9==0 && r10==1 ) || (r7==0 && r8==0 && r9==1 && r10==0 ) || (r7==0 && r8==1 && r9==0 && r10==0 ) || (r7==1 && r8==0 && r9==0 && r10==0 ) || (r7==1 && r8==1 && r9==0 && r10==0 ))
				{
					//RDok[3]++;
					my_cof[i][PLANE_Y][x][y+3] = 0; my_cof[i][PLANE_Y][x+1][y+2] = 0; my_cof[i][PLANE_Y][x+2][y+1] = 0; my_cof[i][PLANE_Y][x+3][y] = 0;

					if(r7==0 && r8==0 && r9==0 && r10==1 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 1; num_sec+=3;}

					if(r7==0 && r8==0 && r9==1 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 0; num_sec+=4;}

					if(r7==0 && r8==1 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; num_sec+=4;}

					if(r7==1 && r8==0 && r9==0 && r10==0 )
					{sec[k++] = 0; sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}
#if 1//PS
					if(r7==1 && r8==1 && r9==0 && r10==0 )
					{sec[k++] = 1; sec[k++] = 1; sec[k++] = 0; sec[k++] = 1; num_sec+=4;}
#endif
					for(numk=0;numk<k;numk++)
						fprintf(fp_sec,"%d",sec[numk]);
					return numcoeff;
				}
			}
			//#endif
#endif


			//#endif
			{
				int start, end;

				if(rdh_encode == 1) { start = 0; end= 4;}
				else if (rdh_encode == 2) { start = 3; end= -1;}

				j= start;
				while(j!=end)
				{
					int dh;
					int x1,y1,x2,y2;

					x1 = x + j;
					x2 = x + ((j+1)%4);
					y1 = y + (3-j);
					y2 = y + (3-((j+1)%4));


					dh = (int)ceil(((double)(my_cof[i][PLANE_Y][x1][y1] - (double)my_cof[i][PLANE_Y][x2][y2])/(double)2));

					//#if RDH_ENCOD
					if(rdh_encode == 1) {
						if(dh == 0)
						{
							char secbit;

							if(randnum[rand_count] == '1')
								secbit = 1;
							else
								secbit = 0;

							my_cof[i][PLANE_Y][x1][y1] = my_cof[i][PLANE_Y][x1][y1] + secbit;
							my_cof[i][PLANE_Y][x2][y2] = my_cof[i][PLANE_Y][x2][y2] - secbit;

							rand_count++;
							rdh_count++;
						}
						else if(dh > 0)
						{
							my_cof[i][PLANE_Y][x1][y1] = my_cof[i][PLANE_Y][x1][y1] + 1;
							my_cof[i][PLANE_Y][x2][y2] = my_cof[i][PLANE_Y][x2][y2] - 1;
						}
					}
					//#else
					else if(rdh_encode == 2) {
						if(dh == 0)
						{
							sec[k] = 0;
							//fprintf(fp_sec,"%d",sec);
							num_sec++;

							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];
							k++;
						}
						else if(dh == 1)
						{
							sec[k] = 1;
							//fprintf(fp_sec,"%d",sec);
							num_sec++;

							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - sec[k];
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + sec[k];

							my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
							my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
							k++;
						}
						else if(dh > 1)
						{
							my_coef[0] = my_cof[i][PLANE_Y][x1][y1] - 1;
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2] + 1;

							my_cof[i][PLANE_Y][x1][y1] = my_coef[0];
							my_cof[i][PLANE_Y][x2][y2] = my_coef[1];
						}
						else if(dh < 0)
						{
							my_coef[0] = my_cof[i][PLANE_Y][x1][y1];
							my_coef[1] = my_cof[i][PLANE_Y][x2][y2];
						}
					}
					//#endif
					if(rdh_encode == 1) j++;
					else if (rdh_encode == 2) j--;
				}
			}



			//#if !RDH_ENCOD
			if (rdh_encode == 2) {
				for(numk=k-1;numk>=0;numk--)
					fprintf(fp_sec,"%d",sec[numk]);
			}
			//#endif
		}
	}		
	return numcoeff;
}



int bcm_rdh_4x4(int****my_cof, int block_x4, int block_y4, int mb_num)  //BOUCHAMA
{
	int numcoeff=0;
	int x,y,i,j;

	int my_coef[2];
	char sec[3];
	int k=0;
	int xx;
	int numk;
	int sec_line;

	i = mb_num;
	x = block_x4;
	y = block_y4;

	if(my_cof[i][PLANE_Y][x+1][y+3]==0 && my_cof[i][PLANE_Y][x+2][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+1]==0 && my_cof[i][PLANE_Y][x+2][y+3]==0 && my_cof[i][PLANE_Y][x+3][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+3]==0)
	{

		{
			//#if RDH_ENCOD
			if(rdh_encode == 1) {
				if(my_cof[i][PLANE_Y][x+3][y+0] == 0 && my_cof[i][PLANE_Y][x+1][y+2] == 0 && my_cof[i][PLANE_Y][x+2][y+1]== 0 && my_cof[i][PLANE_Y][x+0][y+3]== 0)
				{
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '0') ;
					else if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '1')
						my_cof[i][PLANE_Y][x+1][y+2] = 1;
					else if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '0')
						my_cof[i][PLANE_Y][x+2][y+1] = 1;																
					else if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '1')
						my_cof[i][PLANE_Y][x+3][y+0] = 1;
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '0')
						my_cof[i][PLANE_Y][x+0][y+3] = -1; 
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '1')
						my_cof[i][PLANE_Y][x+1][y+2] = -1;
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '0')
						my_cof[i][PLANE_Y][x+2][y+1] = -1;
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '1')
						my_cof[i][PLANE_Y][x+3][y+0] = -1;


					rdh_count += 3;
					rand_count += 3;
				}
				else
				{
					for(j=0;j<4;j++)
					{
						if( (j != 0) && (my_cof[i][PLANE_Y][x+j][y+(3-j)] > 0) )
						{
							my_cof[i][PLANE_Y][x+j][y+(3-j)] +=1;
						}
						else if( my_cof[i][PLANE_Y][x+j][y+(3-j)] < 0 )
						{
							my_cof[i][PLANE_Y][x+j][y+(3-j)] -=1;

						}
					}
				}
			}
			//#else
			else if(rdh_encode == 2) {
				if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 0; sec[1] = 0; sec[2] = 0; 
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3;sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==1 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 0; sec[1] = 0; sec[2] = 1;  my_cof[i][PLANE_Y][x+1][y+2]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3;sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==1 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 0; sec[1] = 1; sec[2] = 0; my_cof[i][PLANE_Y][x+2][y+1]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==1)
				{	sec[0] = 0; sec[1] = 1; sec[2] = 1; my_cof[i][PLANE_Y][x+3][y]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==-1 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 0; sec[2] = 0; my_cof[i][PLANE_Y][x][y+3]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==-1 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 0; sec[2] = 1; my_cof[i][PLANE_Y][x+1][y+2]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==-1 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 1; sec[2] = 0; my_cof[i][PLANE_Y][x+2][y+1]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==-1)
				{	sec[0] = 1; sec[1] = 1; sec[2] = 1; my_cof[i][PLANE_Y][x+3][y]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else
				{
					for(xx=0;xx<4;xx++)
					{
						if((xx!=0) && my_cof[i][PLANE_Y][x+xx][y+(3-xx)] > 0) my_cof[i][PLANE_Y][x+xx][y+(3-xx)]-=1;
						else if(my_cof[i][PLANE_Y][x+xx][y+(3-xx)] < 0) my_cof[i][PLANE_Y][x+xx][y+(3-xx)]+=1;
					}
					sec_line=0;
				}
			}
			//#endif


		}
#if !RDH_ENCOD
		/*if(sec_line == 1)
		{
		for(numk=0;numk<3;numk++)
		fprintf(fp_sec,"%d",sec[numk]);
		}*/
#endif
	}
	return numcoeff;

}



int shd_rdh_4x4(int****my_cof, int block_x4, int block_y4, int mb_num)  //Shaid
{
	int numcoeff=0;
	int x,y,i,j,l;
	char sec;
	int k=0;



	i = mb_num;
	x = block_x4;
	y = block_y4;

	if(my_cof[i][PLANE_Y][x+1][y+0]==0 && my_cof[i][PLANE_Y][x+1][y+1]==0 && my_cof[i][PLANE_Y][x+2][y+0]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0)
		return 0;
	else
	{

		{
			//#if RDH_ENCOD
			if(rdh_encode == 1) {
				for(j=1; j<3 ; j++)
				{
					for(l=0; l<2 ; l++)
					{
						if(my_cof[i][PLANE_Y][x+j][y+l] == 0 )
						{
							sec = randnum[rand_count++];
							rdh_count++;
							my_cof[i][PLANE_Y][x+j][y+l] += atoi(&sec);
						}
						else if(my_cof[i][PLANE_Y][x+j][y+l] > 0 )
							my_cof[i][PLANE_Y][x+j][y+l] +=1;
					}
				}
			}
			//#else
			else if(rdh_encode == 2) {
				for(j=1; j<3 ; j++)
				{
					for(l=0; l<2 ; l++)
					{
						if(my_cof[i][PLANE_Y][x+j][y+l] == 0 )
						{
							sec = 0; num_sec++;
							fprintf(fp_sec,"%d",sec);
						}
						else if(my_cof[i][PLANE_Y][x+j][y+l] == 1 )
						{
							sec = 1; num_sec++;
							my_cof[i][PLANE_Y][x+j][y+l] -=sec;
							fprintf(fp_sec,"%d",sec);
						}
						else if(my_cof[i][PLANE_Y][x+j][y+l] > 1 )
						{
							my_cof[i][PLANE_Y][x+j][y+l] -=1;
						}
					}
				}
			}
			//#endif
		}

	}
	return numcoeff;

}



int kci_rdh_4x4(int****my_cof, int block_x4, int block_y4, int mb_num)  //KCI
{
	int numcoeff=0;
	int x,y,i,j;

	int my_coef[2];
	char sec[4];
	int k=0;
	int xx;
	int numk;
	int sec_line;

	i = mb_num;
	x = block_x4;
	y = block_y4;

	if(my_cof[i][PLANE_Y][x+1][y+3]==0 && my_cof[i][PLANE_Y][x+2][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+1]==0 && my_cof[i][PLANE_Y][x+2][y+3]==0 && my_cof[i][PLANE_Y][x+3][y+2]==0 && my_cof[i][PLANE_Y][x+3][y+3]==0)
	{

		{
			//#if RDH_ENCOD
			if(rdh_encode == 1) {
				if(my_cof[i][PLANE_Y][x+3][y+0] == 0 && my_cof[i][PLANE_Y][x+1][y+2] == 0 && my_cof[i][PLANE_Y][x+2][y+1]== 0 && my_cof[i][PLANE_Y][x+0][y+3]== 0)
				{
					if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '0') 
					{	rand_count += 3; rdh_count += 3;}
					else if(randnum[rand_count] == '0' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '1')
					{	my_cof[i][PLANE_Y][x+1][y+2] = 1;rand_count += 3;rdh_count += 3;}
					else if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '0')
					{	my_cof[i][PLANE_Y][x+2][y+1] = 1;	rand_count += 3;rdh_count += 3;	}														
					else if(randnum[rand_count] == '0' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '1')
					{	my_cof[i][PLANE_Y][x+3][y+0] = 1;rand_count += 3;rdh_count += 3;}
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '0')
					{	my_cof[i][PLANE_Y][x+0][y+3] = -1; rand_count += 3;rdh_count += 3;}
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '0' && randnum[rand_count+2] == '1')
					{	my_cof[i][PLANE_Y][x+1][y+2] = -1;rand_count += 3;rdh_count += 3;}
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '0')
					{	my_cof[i][PLANE_Y][x+2][y+1] = -1;rand_count += 3;rdh_count += 3;}
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '1' && randnum[rand_count+3] == '0')
					{	my_cof[i][PLANE_Y][x+3][y+0] = -1;rand_count += 4;rdh_count += 4;}
					else if(randnum[rand_count] == '1' && randnum[rand_count+1] == '1' && randnum[rand_count+2] == '1' && randnum[rand_count+3] == '1')
					{	my_cof[i][PLANE_Y][x+0][y+3] = 1;rand_count += 4;rdh_count += 4;}


					//rdh_count += 3;
					//rand_count += 3;
				}
				else
				{
					for(j=0;j<4;j++)
					{
						if( my_cof[i][PLANE_Y][x+j][y+(3-j)] > 0 )
						{
							my_cof[i][PLANE_Y][x+j][y+(3-j)] +=1;
						}
						else if( my_cof[i][PLANE_Y][x+j][y+(3-j)] < 0 )
						{
							my_cof[i][PLANE_Y][x+j][y+(3-j)] -=1;

						}
					}
				}
			}
			//#else
			else if(rdh_encode == 2) {
				if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 0; sec[1] = 0; sec[2] = 0; 
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3;sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==1 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 0; sec[1] = 0; sec[2] = 1;  my_cof[i][PLANE_Y][x+1][y+2]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3;sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==1 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 0; sec[1] = 1; sec[2] = 0; my_cof[i][PLANE_Y][x+2][y+1]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==1)
				{	sec[0] = 0; sec[1] = 1; sec[2] = 1; my_cof[i][PLANE_Y][x+3][y]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==-1 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 0; sec[2] = 0; my_cof[i][PLANE_Y][x][y+3]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==-1 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 0; sec[2] = 1; my_cof[i][PLANE_Y][x+1][y+2]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==-1 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 1; sec[2] = 0; my_cof[i][PLANE_Y][x+2][y+1]=0;
				fprintf(fp_sec,"%d%d%d",sec[0],sec[1],sec[2]); num_sec+=3; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==0 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==-1)
				{	sec[0] = 1; sec[1] = 1; sec[2] = 1; sec[3] = 0; my_cof[i][PLANE_Y][x+3][y]=0;
				fprintf(fp_sec,"%d%d%d%d",sec[0],sec[1],sec[2],sec[3]); num_sec+=4; sec_line=1;}
				else if(my_cof[i][PLANE_Y][x][y+3]==1 && my_cof[i][PLANE_Y][x+1][y+2]==0 && my_cof[i][PLANE_Y][x+2][y+1]==0 &&my_cof[i][PLANE_Y][x+3][y]==0)
				{	sec[0] = 1; sec[1] = 1; sec[2] = 1; sec[3] = 1; my_cof[i][PLANE_Y][x][y+3]=0;
				fprintf(fp_sec,"%d%d%d%d",sec[0],sec[1],sec[2],sec[3]); num_sec+=4; sec_line=1;}
				else
				{
					for(xx=0;xx<4;xx++)
					{
						if( my_cof[i][PLANE_Y][x+xx][y+(3-xx)] > 0) my_cof[i][PLANE_Y][x+xx][y+(3-xx)]-=1;
						else if(my_cof[i][PLANE_Y][x+xx][y+(3-xx)] < 0) my_cof[i][PLANE_Y][x+xx][y+(3-xx)]+=1;
					}
					sec_line=0;
				}
			}
			//#endif


		}
#if !RDH_ENCOD
		/*if(sec_line == 1)
		{
		for(numk=0;numk<3;numk++)
		fprintf(fp_sec,"%d",sec[numk]);
		}*/
#endif
	}
	return numcoeff;

}


int copy_rdh_cof4x4(int ***rdh_cof4x4, int num4x4, int *levarr, int *runarr, int numcoeff)
{
	rdh_cof4x4[num4x4];

	return 0;
}

/*!
************************************************************************
* \brief
*    rewrite coefficients (run/level) of 4x4 blocks in a MB
*    from the NAL (CABAC Mode)
************************************************************************
*/
void rewrite_comp_coeff_4x4_CAVLC (Macroblock *currMB, ColorPlane pl, int (*InvLevelScale4x4)[4], int qp_per, int cbp, byte **nzcoeff, int flag)
{
  int block_y, block_x, b8;
  int i, j;
//  int i0, j0;
  int levarr[16] = {0}, runarr[16] = {0};
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;

  int cur_context; 
  int block_y4, block_x4;

  int rw_levarr[17] = {0}, rw_runarr[17] = {0};

  int num4x4=0;

//  EncMacroblock* EncurrMB=NULL;

  if (IS_I16MB(currMB))
  {
	  if (pl == PLANE_Y)
		  cur_context = LUMA_INTRA16x16AC;
	  else if (pl == PLANE_U)
		  cur_context = CB_INTRA16x16AC;
	  else
		  cur_context = CR_INTRA16x16AC;
  }
  else
  {
	  if (pl == PLANE_Y)
		  cur_context = LUMA;
	  else if (pl == PLANE_U)
		  cur_context = CB;
	  else
		  cur_context = CR;
  }

  //확인용 
  /*if(currSlice->current_mb_nr == 272)
	  printf("\n");*/

  for (block_y = 0; block_y < 4; block_y += 2) /* all modes */
  {
	  block_y4 = block_y << 2;
	  for (block_x = 0; block_x < 4; block_x += 2)
	  {
		  block_x4 = block_x << 2;
		  b8 = (block_y + (block_x >> 1));

		  {
			  if(new_cbp & (1<<b8))//rdh한 값(rdh_cof4x4)을 levarr,와 runarr로 만들기
			  {
				  for (j = block_y4; j < block_y4 + 8; j += BLOCK_SIZE)
				  {
					  for (i = block_x4; i < block_x4 + 8; i += BLOCK_SIZE)
					  {
						  //확인용
						  /*if(rw_bs->byte_pos == 23911)
						  {
							  printf("\n");
						  }*/
						  RewirteCopyLevRun(rdh_cof4x4[num4x4][0], rdh_cof4x4[num4x4][1], rw_levarr, rw_runarr, 16);
						  RewriteCoeff4x4_CAVLC_normal(currMB, cur_context, i >> 2, j >> 2, rw_levarr, rw_runarr, rw_bs);
						  num4x4++;
					  }
				  }
			  }
			  else
			  {
				  rdh_nz_coeff[currSlice->current_mb_nr][block_y    ][block_x    ] = 0;
				  rdh_nz_coeff[currSlice->current_mb_nr][block_y    ][block_x + 1] = 0;
				  rdh_nz_coeff[currSlice->current_mb_nr][block_y + 1][block_x    ] = 0;
				  rdh_nz_coeff[currSlice->current_mb_nr][block_y + 1][block_x + 1] = 0;
			  }
		  }
	  }
  }


  //if (pl == PLANE_Y)
  //		 writeCoeff16x16_CAVLC (EncurrMB, PLANE_Y, currMB, levarr, runarr);
  ////else
  // //copy
#if RDH_DE_EXIS
  {
	  if(currSlice->slice_type == I_SLICE)
	  {
		  if(!(currMB->luma_transform_size_8x8_flag))
		  {
			  if(pl == PLANE_Y)
			  {
				  rdh_De4x4_EXIS_CAVLE_restore(currMB, my_cof, InvLevelScale4x4, qp_per);
				  //rdh_De4x4_EXIS_CABAC(currMB, my_cof, InvLevelScale4x4, qp_per);
			  }
		  }
	  }
	  else if(currSlice->slice_type == P_SLICE)
	  {
		  if(!(currMB->luma_transform_size_8x8_flag))
		  {
			  if(pl == PLANE_Y)
			  {
				  rdh_De4x4_EXIS_CAVLE_restore(currMB, my_cof, InvLevelScale4x4, qp_per);	  
			  }
		  }
	  }
  }
#endif 

#if RDH_DE_Double
  {
	  if(currSlice->slice_type == I_SLICE)
	  {
		  if(!(currMB->luma_transform_size_8x8_flag))
		  {
			  if(pl == PLANE_Y)
			  {
				  rdh_De4x4_Double_CAVLE(currMB, my_cof, InvLevelScale4x4, qp_per);	  
			  }
		  }
	  }
	  else if(currSlice->slice_type == P_SLICE)
	  {
		  if(!(currMB->luma_transform_size_8x8_flag))
		  {
			  if(pl == PLANE_Y)
			  {
				  rdh_De4x4_Double_CAVLE_Pframe(currMB, my_cof, InvLevelScale4x4, qp_per);	  
			  }
		  }
	  }
  }
#endif 

#if RDH_DE_PredHS || RDH_DE_PDHS_mod
  {
	  if(currSlice->slice_type == I_SLICE)
	  {
		  if(!(currMB->luma_transform_size_8x8_flag))
		  {
			  if(pl == PLANE_Y)
			  {
				  rdh_De4x4_PredHS_CAVLE_restore(currMB, my_cof, InvLevelScale4x4, qp_per);	  
			  }
		  }
	  }
#if 0
	  else if(currSlice->slice_type == P_SLICE)
	  {
		  if(!(currMB->luma_transform_size_8x8_flag))
		  {
			  if(pl == PLANE_Y)
			  {
				  rdh_De4x4_PredHS_CAVLE_restore(currMB, my_cof, InvLevelScale4x4, qp_per);	  
			  }
		  }
	  }
#endif
  }
#endif 
}

int rw_bit2buf(byte buffer[],int totbitoffset, int numbits)
{
	int bitoffset  = 7 - (totbitoffset & 0x07); // bit from start of byte
    int byteoffset = (totbitoffset >> 3); // byte from start of buffer
    int bitcounter = numbits;
    byte *curbyte  = &(buffer[byteoffset]);
	byte *byte_buf  = &rw_bs->byte_buf;
	int *bits_to_go = &rw_bs->bits_to_go;
//	int i;
	int bitcount= 0;

	 while (numbits--)
    {
      *byte_buf <<=1;    
      *byte_buf |= ((*curbyte)>> (bitoffset--)) & 0x01;    
      if (bitoffset == -1 ) 
      { //Move onto next byte to get all of numbits
        curbyte++;
        bitoffset = 7;
      }
	  if ((--(*bits_to_go)) == 0)
      {
        *bits_to_go = 8;      
        rw_bs->streamBuffer[rw_bs->byte_pos++] = *byte_buf;
        *byte_buf = 0;      
      }  
	  bitcount++;
    }

    return bitcounter;           // return absolute offset in bit from start of frame
}


void  rw_buf2buf(RW_Bitstream *bsBuff)
{
	int numbits = ((bsBuff->stored_byte_pos) << 3) + (8 - bsBuff->stored_bits_to_go);
	int bitoffset  = 7;//7 - (totbitoffset & 0x07); // bit from start of byte
    int byteoffset = 0;//(totbitoffset >> 3); // byte from start of buffer
    int bitcounter = numbits;
	byte *curbyte  = &(bsBuff->streamBuffer[byteoffset]);
	byte *byte_buf  = &rw_bs->byte_buf;
	int *bits_to_go = &rw_bs->bits_to_go;
//	int i;
	int bitcount= 0;

	 while (numbits--)
    {
      *byte_buf <<=1;    
      *byte_buf |= ((*curbyte)>> (bitoffset--)) & 0x01;    
      if (bitoffset == -1 ) 
      { //Move onto next byte to get all of numbits
        curbyte++;
        bitoffset = 7;
      }
	  if ((--(*bits_to_go)) == 0)
      {
        *bits_to_go = 8;      
        rw_bs->streamBuffer[rw_bs->byte_pos++] = *byte_buf;
        *byte_buf = 0;      
      }  
	  bitcount++;
    }

    return;           
}

