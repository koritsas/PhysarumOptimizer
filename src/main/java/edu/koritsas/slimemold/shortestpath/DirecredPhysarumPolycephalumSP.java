package edu.koritsas.slimemold.shortestpath;

import edu.koritsas.slimemold.AbstractDirectedPhysarumPolycephalum;
import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.*;
import org.geotools.graph.structure.basic.BasicDirectedNode;

import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 12/10/2016.
 */
public abstract class DirecredPhysarumPolycephalumSP extends AbstractDirectedPhysarumPolycephalum {
    public DirecredPhysarumPolycephalumSP(DirectedGraph graph, Node sourceNode, Node sinkNode, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, Io, γ, numberOfIterations);
    }
}
