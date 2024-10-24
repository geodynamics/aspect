/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Classes for creation and output of image files
 *
 *  Copyright (C) 2003-2015 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <fstream>
#include <sstream>
#include "ls_image.h"
#include "bstream.h"
#include "font_data.inc"
#include "string_utils.h"
#include "share_utils.h"

using namespace std;

const MP_Font medium_bold_font = { 0, 128, 7, 13, medium_bold_font_data };
const MP_Font giant_font = { 0, 128, 9, 15, giant_font_data };

void Palette::setPredefined (int num)
  {
  clear();
  switch(num)
    {
    case 0:
      add(0,Colour(0,0,0));
      add(1,Colour(1,1,1));
      break;
    case 1:
      add(0,Colour(0,0,0));
      add(0.4f,Colour(0,0,0.5f));
      add(0.75f,Colour(0,0.6f,1));
      add(1,Colour(1,1,1));
      break;
    case 4:
      add(0,Colour(0,0,.5f));
      add(0.15f,Colour(0,0,1));
      add(0.4f,Colour(0,1,1));
      add(0.7f,Colour(1,1,0));
      add(0.9f,Colour(1,.33f,0));
      add(1,Colour(.5f,0,0));
      break;
    case 10: // Planck colours 1
      addb(  0,  0,  0,255);
      addb( 42,  0,112,255);
      addb( 85,  0,221,255);
      addb(127,255,237,217);
      addb(170,255,180,  0);
      addb(212,255, 75,  0);
      addb(255,100,  0,  0);
      break;
    case 11: // Planck colours 2
      addb(  0,  0,  0,255);
      addb( 13, 10, 20,255);
      addb( 26, 30,184,255);
      addb( 38, 80,235,255);
      addb( 52,191,239,250);
      addb( 65,228,240,245);
      addb( 76,241,241,212);
      addb( 77,241,241,212);
      addb( 88,245,240,175);
      addb(101,248,235,130);
      addb(114,250,204, 38);
      addb(127,243,153, 13);
      addb(140,204, 77,  0);
      addb(153,165, 32, 32);
      addb(166,114,  0, 32);
      addb(179,128,128,153);
      addb(192,179,179,204);
      addb(205,204,204,230);
      addb(218,230,230,242);
      addb(231,242,242,250);
      addb(255,252,252,255);
      break;
    case 100: // Comptel
      addb(  0,  0,  0,  0);
      addb( 36, 64,  0,128);
      addb( 71,  0,  0,255);
      addb(117,  0,255,255);
      addb(152,  0,255,  0);
      addb(188,255,255,  0);
      addb(223,255,  0,  0);
      addb(255,255,255,255);
      break;
    default:
      planck_fail("Palette #"+dataToString(num)+" not yet supported.");
    }
  }

void LS_Image::write_char (int xpos, int ypos, const Colour &col, char c,
  int scale)
  {
  planck_assert ((c>=font.offset) && (c<font.offset+font.num_chars),
    "write_char: character out of range");
  for (int i=0; i<font.xpix; ++i)
    for (int j=0; j<font.ypix; ++j)
      {
      int ofs = (c-font.offset)*font.xpix*font.ypix + j*font.xpix + i;
      if (font.data[ofs]!=' ')
        for (int m=0; m<scale; ++m)
          for (int n=0; n<scale; ++n)
            put_pixel(xpos+scale*i+m,ypos+scale*j+n,col);
      }
  }

LS_Image::LS_Image ()
  : font(medium_bold_font) {}

LS_Image::LS_Image (int xres, int yres)
  : font(medium_bold_font), pixel(xres,yres,Colour(0,0,0)) {}

void LS_Image::annotate (int xpos, int ypos, const Colour &col,
  const string &text, int scale)
  {
  for (tsize m=0; m<text.length(); ++m)
    write_char(xpos+m*scale*font.xpix, ypos, col, text[m],scale);
  }

void LS_Image::annotate_centered (int xpos, int ypos, const Colour &col,
  const string &text, int scale)
  {
  xpos-=(scale*text.length()*font.xpix)/2;
  ypos-=scale*font.ypix/2;
  annotate (xpos,ypos,col,text,scale);
  }

void LS_Image::set_font (const MP_Font &fnt)
  { font = fnt; }

void LS_Image::write_TGA (const string &file) const
  {
  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file '" + file + "'");

  tsize xres=pixel.size1(), yres=pixel.size2();

  bostream bo(out);
  const uint8 header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    uint8(xres%256), uint8(xres/256), uint8(yres%256), uint8(yres/256), 24, 32};

  bo.put (header, 18);
  vector <uint8> line(3*xres);
  for (tsize j=0; j<yres; ++j)
    {
    for (tsize i=0; i<xres; ++i)
      {
      line[i*3  ]=pixel[i][j].b;
      line[i*3+1]=pixel[i][j].g;
      line[i*3+2]=pixel[i][j].r;
      }
    bo.put(line.data(),3*xres);
    }
  planck_assert(out,"error writing output file '" + file + "'");
  }

namespace {

void write_equal_range (const arr<Colour8> &px, tsize begin, tsize end,
  vector<uint8> &buf)
  {
  chunkMaker cm (end-begin,128);
  uint64 cbeg, csz;
  while (cm.getNext(cbeg,csz))
    {
    buf.push_back(uint8(csz-1+128));
    buf.push_back(px[begin].b);
    buf.push_back(px[begin].g);
    buf.push_back(px[begin].r);
    }
  }
void write_unequal_range (const arr<Colour8> &px, tsize begin, tsize end,
  vector<uint8> &buf)
  {
  chunkMaker cm (end-begin,128);
  uint64 cbeg, csz;
  while (cm.getNext(cbeg,csz))
    {
    buf.push_back(uint8(csz-1));
    for (tsize cnt=begin+cbeg; cnt< begin+cbeg+csz; ++cnt)
      {
      buf.push_back(px[cnt].b);
      buf.push_back(px[cnt].g);
      buf.push_back(px[cnt].r);
      }
    }
  }

} // unnamed namespace

void LS_Image::write_TGA_rle(const string &file) const
  {
  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file '" + file + "'");

  tsize xres=pixel.size1(), yres=pixel.size2();

  bostream bo(out);
  const uint8 header[18] = { 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    uint8(xres%256), uint8(xres/256), uint8(yres%256), uint8(yres/256), 24, 32};

  bo.put(header,18);
  vector<uint8> buf;
  arr<Colour8> px(xres);
  for (tsize y=0; y<yres; ++y)
    {
    buf.clear();
    for (tsize x=0; x<xres; ++x)
      px[x] = pixel[x][y];
    tsize xstart=0;
    while (xstart<xres)
      {
      if (xstart==xres-1)
        {
        write_unequal_range (px,xstart,xstart+1,buf);
        xstart=xres;
        }
      else
        {
        if (px[xstart+1]==px[xstart]) // range of equal pixels
          {
          tsize xend=xstart+2;
          while ((xend<xres) && (px[xend]==px[xstart])) ++xend;
          write_equal_range (px,xstart,xend,buf);
          xstart=xend;
          }
        else
          {
          tsize xend=xstart+2;
          while ((xend<xres) && (px[xend]!=px[xend-1])) ++xend;
          write_unequal_range (px,xstart,xend,buf);
          xstart=xend;
          }
        }
      }
    bo.put(buf.data(), buf.size());
    }
  planck_assert(out,"error writing output file '" + file + "'");
  }

void LS_Image::write_PPM (const string &file) const
  {
  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file '" + file + "'");

  tsize xres=pixel.size1(), yres=pixel.size2();

  bostream bo(out);

  ostringstream header;
  header << "P6" << endl << xres << endl << yres << endl << 255 << endl;
  string hdrdata = header.str();
  bo.put(hdrdata.c_str(),hdrdata.size());

  vector <uint8> line(3*xres);
  for (tsize j=0; j<yres; ++j)
    {
    for (tsize i=0; i<xres; ++i)
      {
      line[i*3  ]=pixel[i][j].r;
      line[i*3+1]=pixel[i][j].g;
      line[i*3+2]=pixel[i][j].b;
      }
    bo.put(line.data(),3*xres);
    }
  planck_assert(out,"error writing output file '" + file + "'");
  }
