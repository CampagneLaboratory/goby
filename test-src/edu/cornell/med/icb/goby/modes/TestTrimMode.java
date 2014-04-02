/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.modes;

import com.google.protobuf.ByteString;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.lang.MutableString;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 6/7/11
 *         Time: 3:25 PM
 */
public class TestTrimMode {
    @Test
    public void testContains() {

        TrimMode trimmer=new TrimMode();
        MutableString original=new MutableString("12>>ACT1234567");
        byte[] bytes=new byte[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
        ByteArrayList list=new ByteArrayList();
        MutableString[] adapters= {
                new MutableString("ACT")
        };
        MutableString result = trimmer.contains(10, original, ByteString.copyFrom(bytes), list, adapters);
        assertEquals(new MutableString("1234567"),result);
        assertEquals(ByteArrayList.wrap(new byte[]{7,8,9,10,11,12,13,14}), list);


    }


    @Test
    public void testLeft() {

        TrimMode trimmer=new TrimMode();
        MutableString original=new MutableString("ACT1234567");
        byte[] bytes=new byte[]{0,1,2,3,4,5,6,7,8,9,10};
        ByteArrayList list=new ByteArrayList();
        MutableString[] adapters= {
                new MutableString("ACT")
        };
        MutableString result = trimmer.trimLeft(10, original, ByteString.copyFrom(bytes), list, adapters);
        assertEquals(new MutableString("1234567"),result);
        assertEquals(ByteArrayList.wrap(new byte[]{3,4,5,6,7,8,9,10}), list);


    }
     @Test
    public void testRight() {

        TrimMode trimmer=new TrimMode();
        MutableString original=new MutableString("ACT1234567ACT");
        byte[] bytes=new byte[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13};
        ByteArrayList list=new ByteArrayList();
        MutableString[] adapters= {
                new MutableString("ACT")
        };
        MutableString result = trimmer.trimRight(original.length(), original, ByteString.copyFrom(bytes), list, adapters);
        assertEquals(new MutableString("ACT1234567"),result);
        assertEquals(ByteArrayList.wrap(new byte[]{0,1,2,3,4,5,6,7,8,9}), list);


    }

     @Test
    public void testMinTrimLengthLeft() {

        TrimMode trimmer=new TrimMode();
         trimmer.setMinLengthLeft(4);
        MutableString original=new MutableString("ACT1234567");
        byte[] bytes=new byte[]{0,1,2,3,4,5,6,7,8,9,10};
        ByteArrayList list=new ByteArrayList();
        MutableString[] adapters= {
                new MutableString("ACT")
        };
        MutableString result = trimmer.trimLeft(10, original, ByteString.copyFrom(bytes), list, adapters);
        assertEquals(new MutableString("ACT1234567"),result);


    }
    @Test
    public void testMinLengthRight() {

        TrimMode trimmer=new TrimMode();
        trimmer.setMinLengthRight(4);
        MutableString original=new MutableString("ACT1234567ACT");
        byte[] bytes=new byte[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13};
        ByteArrayList list=new ByteArrayList();
        MutableString[] adapters= {
                new MutableString("ACT")
        };
        MutableString result = trimmer.trimRight(original.length(), original, ByteString.copyFrom(bytes), list, adapters);
        assertEquals(new MutableString("ACT1234567ACT"),result);


    }

}
