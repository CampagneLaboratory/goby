package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.util.Collection;
import java.util.Map;

/**
 * @author Fabien Campagne
 *         Date: 1/26/13
 *         Time: 12:37 PM
 */
public class PositionToBasesMap<T> {
    private Int2ObjectAVLTreeMap<T> delegate = new Int2ObjectAVLTreeMap<T>();
    private IntSortedSet sortedKeys = new IntAVLTreeSet();
    private Int2BooleanAVLTreeMap ignoredPositions = new Int2BooleanAVLTreeMap();

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        IntSortedSet sorted = new IntAVLTreeSet();
        sorted.addAll(delegate.keySet());

        builder.append(String.format("key span: [%d-%d]%n", sorted.firstInt(), sorted.lastInt()));
        for (T value : delegate.values()) {
            builder.append(value.toString());
            builder.append("\n");
        }
        return builder.toString();
    }

    public IntSet keySet() {
        return sortedKeys;
    }

    public boolean containsKey(int k) {
        return delegate.containsKey(k);
    }

    public int size() {
        return delegate.size();
    }

    public boolean isEmpty() {
        return delegate.isEmpty();
    }

    public void clear() {
        sortedKeys.clear();
        delegate.clear();
        ignoredPositions.clear();
    }

    public T remove(int k) {
        sortedKeys.remove(k);
        ignoredPositions.remove(k);
        return delegate.remove(k);

    }

    public T get(int ok) {
        return delegate.get(ok);
    }

    public void put(int keyPos, T positionBaseInfos) {
        sortedKeys.add(keyPos);
        delegate.put(keyPos, positionBaseInfos);
    }

    public ObjectSet<Map.Entry<Integer, T>> entrySet() {
        return delegate.entrySet();
    }

    public int firstPosition() {
        return sortedKeys.firstInt();
    }

    public void markIgnoredPosition(int position) {
        ignoredPositions.put(position,true);
    }
    public boolean isIgnoredPosition(int position) {
        return ignoredPositions.get(position);
    }
}
