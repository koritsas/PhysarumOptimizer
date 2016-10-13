package edu.koritsas.slimemold;

import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.shortestpath.DirecredPhysarumPolycephalumSP;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

/**
 * Created by ilias on 13/10/2016.
 */
public class DirectedSlimeSP extends DirecredPhysarumPolycephalumSP {
    public DirectedSlimeSP(DirectedGraph graph, Node source, Node sink, double Io, double γ, int numberOfIterations) {
        super(graph, source, sink, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeCost(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry g = (Geometry) f.getDefaultGeometry();
        return g.getLength();
    }
}
