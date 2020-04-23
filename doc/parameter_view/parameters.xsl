<?xml version="1.0" encoding="UTF-8"?>

<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://www.w3.org/1999/xhtml">

	<xsl:template match="/">
		<html>
			<head>
				<link rel="stylesheet" type="text/css" href="parameters.css"></link>
			</head>
			<body>
				<h2>
					ASPECT input parameters  
				</h2>

				<div>
					<button onclick="expandAll()">Expand all</button> 
					<button onclick="collapseAll()">Collapse all</button> 
					<button onclick="expandAllSubsections()">Expand subsections</button> 
					<button onclick="collapseAllSubsections()">Collapse subsections</button>
					<input type="text" id="search" value=""/>
					<a href="parameters.xml" id="link">link</a>
				</div>
				<ul id="ParameterList">
					<xsl:apply-templates
						select="ParameterHandler/*" />
				</ul>
                                <script src="parameters.js" type="text/javascript"></script>
			</body>
		</html>
	</xsl:template>

	<xsl:template match="ParameterHandler//*">

			<xsl:choose>
				<xsl:when test=".//value">
					<li>
						                        <xsl:choose>
                                <xsl:when test="./value">
						<div class="collapsible parameter mangled">
							set <b> <xsl:value-of select="name()" /> </b> = <em> <xsl:value-of select="./default_value" /> </em> <xsl:value-of select="./pattern_description" />
						</div>
						        </xsl:when>
							<xsl:otherwise>
		                                <div class="collapsible subsection mangled">
							subsection <b> <xsl:value-of select="name()" /> </b>
                                                </div>
					</xsl:otherwise>
				</xsl:choose>
			<div class="content">
					<ul>
						<xsl:apply-templates select="child::*" />
					</ul>
				</div>
		</li>
	</xsl:when>
	<xsl:when test="name() = 'value'"></xsl:when>
	<xsl:when test="name() = 'default_value'"></xsl:when>
	<xsl:when test="name() = 'pattern'"></xsl:when>
	<xsl:when test="name() = 'pattern_description'"></xsl:when>
				<xsl:otherwise>
		<li>
					<xsl:value-of select="." />
		</li>
				</xsl:otherwise>
			</xsl:choose>
	</xsl:template>
</xsl:stylesheet>
