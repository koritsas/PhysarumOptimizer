package edu.koritsas.slimemold;

import com.vividsolutions.jts.geom.Geometry;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

/**
 * Created by ilias on 12/10/2016.
 */
public class DirectedSlimeSp extends DirecredPhysarumPolycephalumSP {
    public DirectedSlimeSp(DirectedGraph graph, Node source, Node sink, double Io, double γ, int numberOfIterations) {
        super(graph, source, sink, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeCost(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry g = (Geometry) f.getDefaultGeometry();
        return g.getLength();
    }
}
