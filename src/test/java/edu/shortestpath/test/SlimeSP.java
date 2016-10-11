package edu.shortestpath.test;

import com.vividsolutions.jts.geom.Geometry;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

/**
 * Created by ilias on 6/10/2016.
 */
public class SlimeSP extends  PhysarumPolycephalumSP{
    public SlimeSP(Graph graph, Node sourceNode, Node sinkNode, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeCost(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry geom = (Geometry) f.getDefaultGeometry();
        return geom.getLength();
    }
}
