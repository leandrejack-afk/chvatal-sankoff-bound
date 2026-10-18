// Reviewer tool: check rows of the new C table trans_h<h>.bin against Lueker's own
// State.transition (unmodified State/DPTable classes) on random states.
// usage: java CheckWithLueker <trans file> <samples> <seed>
import java.io.*;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.util.Random;

public class CheckWithLueker {
    public static void main(String[] a) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(a[0], "r");
        FileChannel ch = raf.getChannel();
        ByteBuffer hb = ByteBuffer.allocate(8).order(ByteOrder.LITTLE_ENDIAN);
        ch.read(hb, 0); hb.flip();
        int h = hb.getInt(); int n = hb.getInt();
        State.h = h;
        int L = 2 * h - 1, NT = 2 * L;
        long koff = 8L + 4L * NT * n;
        long[] keys = new long[n];
        ByteBuffer kb = ByteBuffer.allocate(1 << 24).order(ByteOrder.LITTLE_ENDIAN);
        long pos = koff; int idx = 0;
        while (idx < n) {
            kb.clear();
            int want = (int) Math.min(kb.capacity(), 8L * (n - idx));
            kb.limit(want);
            while (kb.hasRemaining()) { int r = ch.read(kb, pos + kb.position()); if (r < 0) throw new IOException("eof"); }
            kb.flip();
            while (kb.remaining() >= 8) keys[idx++] = kb.getLong();
            pos += want;
        }
        int mask = (1 << h) - 1;
        int samples = Integer.parseInt(a[1]);
        Random rnd = new Random(Long.parseLong(a[2]));
        ByteBuffer rb = ByteBuffer.allocate(4 * NT).order(ByteOrder.LITTLE_ENDIAN);
        long mism = 0, checked = 0, acc = 0;
        for (int t = 0; t < samples; t++) {
            int s = (t == 0) ? 1 : 1 + rnd.nextInt(n);
            rb.clear();
            while (rb.hasRemaining()) ch.read(rb, 8L + 4L * NT * (s - 1) + rb.position());
            rb.flip();
            long k = keys[s - 1];
            State st = new State((int) ((k >> (3 * h)) & mask), (int) ((k >> (2 * h)) & mask),
                                 (int) ((k >> h) & mask), (int) (k & mask));
            for (int c = 0; c < 2; c++) for (int l = -(h - 1); l <= h - 1; l++) {
                int tgt = rb.getInt();
                State ns = st.transition(c, l);
                boolean ok;
                if (ns == null) ok = (tgt == 0);
                else {
                    acc++;
                    if (tgt < 1 || tgt > n) ok = false;
                    else {
                        long kk = keys[tgt - 1];
                        ok = ((kk >> (3 * h)) & mask) == (ns.x & 0xffff) && ((kk >> (2 * h)) & mask) == (ns.y & 0xffff)
                          && ((kk >> h) & mask) == (ns.dx & 0xffff) && (kk & mask) == (ns.dy & 0xffff);
                    }
                }
                if (!ok) { mism++; if (mism < 10) System.out.println("mismatch s=" + s + " c=" + c + " l=" + l); }
                checked++;
            }
        }
        System.out.println("h=" + h + " N=" + n + " sampled states=" + samples + " transitions checked=" + checked
                           + " accepted=" + acc + " mismatches vs Lueker State.transition=" + mism);
    }
}
