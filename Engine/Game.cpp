/****************************************************************************************** 
 *	Chili DirectX Framework Version 16.07.20											  *	
 *	Game.cpp																			  *
 *	Copyright 2016 PlanetChili.net <http://www.planetchili.net>							  *
 *																						  *
 *	This file is part of The Chili DirectX Framework.									  *
 *																						  *
 *	The Chili DirectX Framework is free software: you can redistribute it and/or modify	  *
 *	it under the terms of the GNU General Public License as published by				  *
 *	the Free Software Foundation, either version 3 of the License, or					  *
 *	(at your option) any later version.													  *
 *																						  *
 *	The Chili DirectX Framework is distributed in the hope that it will be useful,		  *
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of						  *
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the						  *
 *	GNU General Public License for more details.										  *
 *																						  *
 *	You should have received a copy of the GNU General Public License					  *
 *	along with The Chili DirectX Framework.  If not, see <http://www.gnu.org/licenses/>.  *
 ******************************************************************************************/
#include "MainWindow.h"
#include "Game.h"
#include <cassert>

Game::Game( MainWindow& wnd )
	:
	wnd( wnd ),
	gfx( wnd )
{
}

void Game::Go()
{
	gfx.BeginFrame();	
	UpdateModel();
	ComposeFrame();
	gfx.EndFrame();
}

void Game::UpdateModel()
{
	if (wnd.mouse.LeftIsPressed())
	{
		Vec2 mousePos(float(wnd.mouse.GetPosX()), float(wnd.mouse.GetPosY()));
		mousePos /= cellDimension;
		density[GetId(int(mousePos.x), int(mousePos.y))] += 0.1f;
		if (density[GetId(int(mousePos.x), int(mousePos.y))] > 1.0f)
		{
			density[GetId(int(mousePos.x), int(mousePos.y))] = 1.0f;
		}
	}
}

int Game::GetId(int i, int j)
{
	return (j+1) * (n+2) + (i+1);
}

void Game::DrawDensity()
{
	//Drawing every cell
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			const int A = int(j * cellDimension);
			const int B = int(i * cellDimension);
			for (int J = A; J < A+cellDimension; J++)
			{
				for (int I = B; I < B+cellDimension; I++)
				{
					const int id = GetId(i, j);
					gfx.PutPixel(I, J, Color( char(255 * density[id]), char(255 * density[id]), char(255 * density[id])));
				}
			}
		}
	}
}

void Game::DrawVelocities(bool separated)
{
	//Given that all velocities are normalized we can use this function to draw the vectors
	const float maxDrawnLength = cellDimension / 2.0f;
	if (separated)
	{
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				const Vec2 c(cellDimension * i + maxDrawnLength, cellDimension * j + maxDrawnLength);
				//Horizontal lines
				gfx.DrawLine(c, c + Vec2(velocity[GetId(i,j)].x, 0.0f) * maxDrawnLength, Colors::Green);
				//Vertical lines
				gfx.DrawLine(c, c + Vec2(0.0f, velocity[GetId(i, j)].y) * maxDrawnLength, Colors::Red);
			}
		}
	}
	else
	{
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				const Vec2 c(cellDimension * i + maxDrawnLength, cellDimension * j + maxDrawnLength);
				gfx.DrawLine(c, c + velocity[GetId(i, j)] * maxDrawnLength, Colors::White);
			}
		}
	}
}

void Game::ComposeFrame()
{
	DrawDensity();
}
